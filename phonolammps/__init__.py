__version__ = '0.4.3'

import numpy as np

from lammps import lammps
from phonopy.file_IO import write_FORCE_CONSTANTS, write_force_constants_to_hdf5
from phonolammps.arrange import get_correct_arrangement
from phonolammps.phonopy_link import obtain_phonon_dispersion_bands, get_phonon
from phonolammps.iofile import get_structure_from_poscar, get_structure_from_lammps, generate_VASP_structure


class Phonolammps:
    def __init__(self,
                 lammps_input_file,
                 supercell_matrix=np.identity(3),
                 primitive_matrix=np.identity(3),
                 displacement_distance=0.01,
                 show_log=False,
                 show_progress=False):
        """
        :param lammps_input_file:  LAMMPS input file name (see example)

        :param displacement_distance: displacement distance in Angstroms
        :return: data_sets phonopy dictionary (without forces) and list of numpy arrays containing
        """

        self._structure = get_structure_from_lammps(lammps_input_file)
        self._lammps_input_file = lammps_input_file

        self._supercell_matrix = supercell_matrix
        self._primitive_matrix = primitive_matrix
        self._displacement_distance = displacement_distance
        self._show_log = show_log
        self._show_progress = show_progress

        self._force_constants = None

    def get_lammps_forces(self, cell_with_disp):
        """
        Calculate the forces of a supercell using lammps
        :param cell_with_disp: supercell from which determine the forces
        :return: numpy array matrix with forces of atoms [Natoms x 3]
        """

        structure = self._structure
        supercell_matrix = self._supercell_matrix
        input_file = self._lammps_input_file

        cmd_list = ['-log', 'none']
        if not self._show_log:
            cmd_list += ['-echo', 'none', '-screen', 'none']

        lmp = lammps(cmdargs=cmd_list)

        supercell_sizes = np.diag(supercell_matrix)

        lmp.file(input_file)
        lmp.command('replicate {} {} {}'.format(*supercell_sizes))
        lmp.command('run 0')

        na = lmp.get_natoms()

        xc = lmp.gather_atoms("x", 1, 3)
        reference = np.array([xc[i] for i in range(na * 3)]).reshape((na, 3))

        template = get_correct_arrangement(reference, structure)
        indexing = np.argsort(template)

        for i in range(na):
            lmp.command('set atom {} x {} y {} z {}'.format(i + 1,
                                                            cell_with_disp[template[i], 0],
                                                            cell_with_disp[template[i], 1],
                                                            cell_with_disp[template[i], 2]))

        lmp.command('run 0')

        forces = lmp.gather_atoms("f", 1, 3)

        forces = np.array([forces[i] for i in range(na * 3)]).reshape((na, 3))[indexing, :]

        lmp.close()

        return forces

    def get_path_using_seek_path(self):

        """
        Obtain the path in reciprocal space to plot the phonon band structure

        :return: dictionary with list of q-points and labels of high symmetry points
        """

        try:
            import seekpath

            structure = self._structure
            cell = structure.get_cell()
            positions = structure.get_scaled_positions()
            numbers = np.unique(structure.get_chemical_symbols(), return_inverse=True)[1]
            structure = (cell, positions, numbers)
            path_data = seekpath.get_path(structure)

            labels = path_data['point_coords']

            band_ranges = []
            for set in path_data['path']:
                band_ranges.append([labels[set[0]], labels[set[1]]])

            return {'ranges': band_ranges,
                    'labels': path_data['path']}
        except ImportError:
            print ('Seekpath not installed. Autopath is deactivated')
            band_ranges = ([[[0.0, 0.0, 0.0], [0.5, 0.0, 0.5]]])
            return {'ranges': band_ranges,
                    'labels': [['GAMMA', '1/2 0 1/2']]}

    def get_force_constants(self):
        """
        calculate the force constants with phonopy using lammps to calculate forces

        :return: ForceConstants type object containing force constants
        """

        if self._force_constants is None:
            phonon = get_phonon(self._structure,
                                setup_forces=False,
                                super_cell_phonon=self._supercell_matrix)

            phonon.get_displacement_dataset()
            phonon.generate_displacements(distance=self._displacement_distance)
            cells_with_disp = phonon.get_supercells_with_displacements()

            cells_with_disp = [cell.get_positions() for cell in cells_with_disp]

            data_sets = phonon.get_displacement_dataset()

            # Get forces from lammps
            for i, cell in enumerate(cells_with_disp):
                if self._show_progress:
                    print ('displacement {} / {}'.format(i+1, len(cells_with_disp)))
                forces = self.get_lammps_forces(cell)
                data_sets['first_atoms'][i]['forces'] = forces

            phonon.set_displacement_dataset(data_sets)
            phonon.produce_force_constants()
            self._force_constants = phonon.get_force_constants()

        return self._force_constants

    def plot_phonon_dispersion_bands(self):
        """
        Plot phonon band structure using seekpath automatic k-path
        """

        import matplotlib.pyplot as plt

        def replace_list(text_string):
            substitutions = {'GAMMA': u'$\Gamma$',
                             }

            for item in substitutions.items():
                text_string = text_string.replace(item[0], item[1])
            return text_string

        force_constants = self.get_force_constants()
        bands_and_labels = self.get_path_using_seek_path()

        _bands = obtain_phonon_dispersion_bands(self._structure,
                                                bands_and_labels['ranges'],
                                                force_constants,
                                                self._supercell_matrix,
                                                primitive_matrix=self._primitive_matrix,
                                                band_resolution=30)

        for i, freq in enumerate(_bands[1]):
            plt.plot(_bands[1][i], _bands[2][i], color='r')

            # plt.axes().get_xaxis().set_visible(False)
        plt.axes().get_xaxis().set_ticks([])

        plt.ylabel('Frequency [THz]')
        plt.xlabel('Wave vector')
        plt.xlim([0, _bands[1][-1][-1]])
        plt.axhline(y=0, color='k', ls='dashed')
        plt.suptitle('Phonon dispersion')

        if 'labels' in bands_and_labels:
            plt.rcParams.update({'mathtext.default': 'regular'})

            labels = bands_and_labels['labels']

            labels_e = []
            x_labels = []
            for i, freq in enumerate(_bands[1]):
                if labels[i][0] == labels[i - 1][1]:
                    labels_e.append(replace_list(labels[i][0]))
                else:
                    labels_e.append(
                        replace_list(labels[i - 1][1]) + '/' + replace_list(labels[i][0]))
                x_labels.append(_bands[1][i][0])
            x_labels.append(_bands[1][-1][-1])
            labels_e.append(replace_list(labels[-1][1]))
            labels_e[0] = replace_list(labels[0][0])

            plt.xticks(x_labels, labels_e, rotation='horizontal')

        plt.show()

    def write_force_constants(self, filename='FORCE_CONSTANTS', hdf5=False):
        """
        Write the force constants in a file in phonopy plain text format

        :param filename: Force constants filename
        :return:
        """

        force_constants = self.get_force_constants()
        if hdf5:
            write_force_constants_to_hdf5(force_constants, filename=filename)
        else:
            write_FORCE_CONSTANTS(force_constants, filename=filename)

    def get_unitcell(self):
        return self._structure

    def get_supercell_matrix(self):
        return self._supercell_matrix

    def get_seekpath_bands(self, band_resolution=30):
        ranges = self.get_path_using_seek_path()['ranges']
        bands =[]
        for q_start, q_end in ranges:
            band = []
            for i in range(band_resolution+1):
                band.append(np.array(q_start) + (np.array(q_end) - np.array(q_start)) / band_resolution * i)
            bands.append(band)

        return bands

    def write_unitcell_POSCAR(self, filename='POSCAR'):
        """
        Write unit cell in VASP POSCAR type file

        :param filename:
        :return:
        """
        poscar_txt = generate_VASP_structure(self._structure)

        with open(filename, mode='w') as f:
            f.write(poscar_txt)
