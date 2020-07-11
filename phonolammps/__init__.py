__version__ = '0.6.0'

import numpy as np

from phonopy.file_IO import write_FORCE_CONSTANTS, write_force_constants_to_hdf5, write_FORCE_SETS
from phonolammps.arrange import get_correct_arrangement, rebuild_connectivity_tinker
from phonolammps.phonopy_link import obtain_phonon_dispersion_bands, get_phonon
from phonolammps.iofile import get_structure_from_poscar, generate_VASP_structure
from phonolammps.iofile import generate_tinker_key_file, generate_tinker_txyz_file, parse_tinker_forces
from phonolammps.iofile import get_structure_from_lammps, get_structure_from_txyz

# define the force unit conversion factors to LAMMPS metal style (eV/Angstrom)
unit_factors = {'real': 4.336410389526464e-2,
                'metal': 1.0,
                'si': 624150636.3094,
                'tinker': 0.043  # kcal/mol to eV
}


class PhonoBase:
    """
    Base class for PhonoLAMMPS
    This class is not designed to be called directly.
    To use it make a subclass and implement the following methods:

    * __init__()
    * get_forces()

    """

    def get_path_using_seek_path(self):

        """ Obtain the path in reciprocal space to plot the phonon band structure

        :return: dictionary with list of q-points and labels of high symmetry points
        """

        try:
            import seekpath

            cell = self._structure.get_cell()
            positions = self._structure.get_scaled_positions()
            numbers = np.unique(self._structure.get_chemical_symbols(), return_inverse=True)[1]

            path_data = seekpath.get_path((cell, positions, numbers))

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

    def get_force_constants(self, include_data_set=False):
        """
        calculate the force constants with phonopy using lammps to calculate forces

        :return: ForceConstants type object containing force constants
        """

        if self._force_constants is None:
            phonon = get_phonon(self._structure,
                                setup_forces=False,
                                super_cell_phonon=self._supercell_matrix,
                                primitive_matrix=self._primitive_matrix,
                                NAC=self._NAC,
                                symmetrize=self._symmetrize)

            phonon.get_displacement_dataset()
            phonon.generate_displacements(distance=self._displacement_distance)
            cells_with_disp = phonon.get_supercells_with_displacements()

            #cells_with_disp = [cell.get_positions() for cell in cells_with_disp]

            data_set = phonon.get_displacement_dataset()

            # Get forces from lammps
            for i, cell in enumerate(cells_with_disp):
                if self._show_progress:
                    print('displacement {} / {}'.format(i+1, len(cells_with_disp)))
                forces = self.get_forces(cell)
                data_set['first_atoms'][i]['forces'] = forces

            phonon.set_displacement_dataset(data_set)
            phonon.produce_force_constants()
            self._force_constants = phonon.get_force_constants()
            self._data_set = data_set

        if include_data_set:
            return [self._force_constants, self._data_set]
        else:
            return self._force_constants

    def plot_phonon_dispersion_bands(self):
        """
        Plot phonon band structure using seekpath automatic k-path
        Warning: The labels may be wrong if the structure is not standarized

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
        """

        force_constants = self.get_force_constants()
        if hdf5:
            write_force_constants_to_hdf5(force_constants, filename=filename)
        else:
            write_FORCE_CONSTANTS(force_constants, filename=filename)

    def write_force_sets(self, filename='FORCE_SETS'):
        """
        Write the force sets in a file in phonopy plain text format

        :param filename: Force sets filename
        """

        data_set = self.get_force_constants(include_data_set=True)[1]

        write_FORCE_SETS(data_set, filename=filename)

    def get_unitcell(self):
        """
        Get unit cell structure

        :return unitcell: unit cell 3x3 matrix (lattice vectors in rows)
        """
        return self._structure

    def get_supercell_matrix(self):
        """
        Get the supercell matrix

        :return supercell: the supercell 3x3 matrix (list of lists)
        """
        return self._supercell_matrix

    def get_primitve_matrix(self):
        return self._primitive_matrix

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

        :param filename: POSCAR file name (Default: POSCAR)
        """
        poscar_txt = generate_VASP_structure(self._structure)

        with open(filename, mode='w') as f:
            f.write(poscar_txt)


################################
#            LAMMPS            #
################################
class Phonolammps(PhonoBase):
    def __init__(self,
                 lammps_input,
                 supercell_matrix=np.identity(3),
                 primitive_matrix=np.identity(3),
                 displacement_distance=0.01,
                 show_log=False,
                 show_progress=False,
                 use_NAC=False,
                 symmetrize=True):
        """
        Main PhonoLAMMPS class

        :param lammps_input: LAMMPS input file name or list of commands
        :param supercell_matrix:  3x3 matrix supercell
        :param primitive cell:  3x3 matrix primitive cell
        :param displacement_distance: displacement distance in Angstroms
        :param show_log: Set true to display lammps log info
        :param show_progress: Set true to display progress of calculation
        """

        # Check if input is file or list of commands
        if type(lammps_input) is str:
            # read from file name
            self._lammps_input_file = lammps_input
            self._lammps_commands_list = open(lammps_input).read().split('\n')
        else:
            # read from commands
            self._lammps_commands_list = lammps_input

        self._structure = get_structure_from_lammps(self._lammps_commands_list)

        self._supercell_matrix = supercell_matrix
        self._primitive_matrix = primitive_matrix
        self._displacement_distance = displacement_distance
        self._show_log = show_log
        self._show_progress = show_progress
        self._symmetrize = symmetrize
        self._NAC = use_NAC

        self._force_constants = None
        self._data_set = None

        self.units = self.get_units(self._lammps_commands_list)

        if not self.units in unit_factors.keys():
            print ('Units style not supported, use: {}'.format(unit_factors.keys()))
            exit()

    def get_units(self, commands_list):
        """
        Get the units label for LAMMPS "units" command from a list of LAMMPS input commands

        :param commands_list: list of LAMMPS input commands (strings)
        :return units: string containing the units
        """
        for line in commands_list:
                if line.startswith('units'):
                    return line.split()[1]
        return 'lj'

    def get_forces(self, cell_with_disp):
        """
        Calculate the forces of a supercell using lammps

        :param cell_with_disp: supercell from which determine the forces
        :return: numpy array matrix with forces of atoms [Natoms x 3]
        """

        import lammps

        supercell_sizes = np.diag(self._supercell_matrix)

        cmd_list = ['-log', 'none']
        if not self._show_log:
            cmd_list += ['-echo', 'none', '-screen', 'none']

        lmp = lammps.lammps(cmdargs=cmd_list)
        lmp.commands_list(self._lammps_commands_list)
        lmp.command('replicate {} {} {}'.format(*supercell_sizes))
        lmp.command('run 0')

        na = lmp.get_natoms()
        xc = lmp.gather_atoms("x", 1, 3)
        reference = np.array([xc[i] for i in range(na * 3)]).reshape((na, 3))

        template = get_correct_arrangement(reference, self._structure, self._supercell_matrix)
        indexing = np.argsort(template)

        coordinates = cell_with_disp.get_positions()

        for i in range(na):
            lmp.command('set atom {} x {} y {} z {}'.format(i + 1,
                                                            coordinates[template[i], 0],
                                                            coordinates[template[i], 1],
                                                            coordinates[template[i], 2]))

        lmp.command('run 0')

        forces = lmp.gather_atoms("f", 1, 3)

        forces = np.array([forces[i] for i in range(na * 3)]).reshape((na, 3))[indexing, :]
        forces = forces * unit_factors[self.units]
        lmp.close()

        return forces


################################
#            TINKER            #
################################
class PhonoTinker(PhonoBase):

    def __init__(self,
                 txyz_input_file,
                 key_input_file,
                 force_field_file,
                 supercell_matrix=np.identity(3),
                 primitive_matrix=np.identity(3),
                 displacement_distance=0.01,
                 show_log=False,
                 show_progress=False,
                 use_NAC=False,
                 symmetrize=True):
        """
        Experimental class to use Tinker to calculate forces, can be used
        as an example to how to expand phonoLAMMPS to other software

        :param txyz_input_file:  TXYZ input file name (see example)
        :param supercell_matrix:  3x3 matrix supercell
        :param primitive cell:  3x3 matrix primitive cell
        :param displacement_distance: displacement distance in Angstroms
        :param show_log: set true to display lammps log info
        :param show_progress: set true to display progress of calculation
        :param use_NAC: set true to use Non-Analytical corrections or not
        :param symmetrize: set true to use symmetrization of the force constants
        """

        self._structure = get_structure_from_txyz(txyz_input_file, key_input_file)
        self._txyz_input_file = txyz_input_file

        self._supercell_matrix = supercell_matrix
        self._primitive_matrix = primitive_matrix
        self._displacement_distance = displacement_distance
        self._show_log = show_log
        self._show_progress = show_progress
        self._symmetrize = symmetrize
        self._NAC = use_NAC

        self._force_constants = None
        self._data_set = None

        self.units = 'tinker'

        self.force_field = force_field_file

        if not self.units in unit_factors.keys():
            print ('Units style not supported, use: {}'.format(unit_factors.keys()))
            exit()

    def get_forces(self, cell_with_disp):
        """
        Calculate the forces of a supercell using tinker
        :param cell_with_disp: supercell (PhonopyAtoms) from which determine the forces
        :return array: numpy array matrix with forces of atoms [Natoms x 3]
        """

        import tempfile
        import subprocess
        import os
        from subprocess import PIPE

        temp_file_name = tempfile.gettempdir() + '/tinker_temp' + '_' + str(os.getpid())

        # temp_file_name = 'test_calc'

        supercell_wd = rebuild_connectivity_tinker(self._structure,
                                                   cell_with_disp,
                                                   self._supercell_matrix)

        tinker_input_file = open(temp_file_name + '.txyz', mode='w')
        tinker_input_file.write(generate_tinker_txyz_file(supercell_wd))

        tinker_key_file = open(temp_file_name + '.key', mode='w')
        tinker_key_file.write(generate_tinker_key_file(supercell_wd))

        tinker_input_file.close()
        tinker_key_file.close()

        tinker_command = './testgrad ' + tinker_input_file.name + \
                         ' ' + self.force_field + ' Y N N ' + ' -k ' + tinker_key_file.name

        tinker_process = subprocess.Popen(tinker_command, stdin=PIPE, stderr=PIPE, stdout=PIPE, shell=True)
        (output, err) = tinker_process.communicate()
        tinker_process.wait()

        if len(err.split()) != 0:
            print(err)
            print('Something wrong in forces calculation!')
            exit()

        # print(output)
        os.unlink(tinker_input_file.name)
        os.unlink(tinker_key_file.name)

        forces = parse_tinker_forces(output) * unit_factors[self.units]

        return forces


if __name__ == '__main__':

    structure = get_structure_from_txyz('structure_wrap_min.txyz', 'structure.key')
    print(structure)
    print(structure.get_connectivity())
    print(generate_VASP_structure(structure))
    print(structure.get_scaled_positions())
    print(structure.get_chemical_symbols())

    phonon = get_phonon(structure,
                        setup_forces=False,
                        super_cell_phonon=[[2, 0, 0], [0, 2, 0], [0, 0, 2]],
                        NAC=False,
                        symmetrize=True)


    phonon.get_displacement_dataset()
    phonon.generate_displacements(distance=0.0001)
    cells_with_disp = phonon.get_supercells_with_displacements()
    print(cells_with_disp[0])
    print(generate_VASP_structure(cells_with_disp[0]))

    supercell_wd = rebuild_connectivity_tinker(structure,
                                               cells_with_disp[0],
                                               phonon.get_supercell_matrix())

    print(generate_tinker_txyz_file(supercell_wd))
    print(generate_tinker_key_file(supercell_wd))

    print(generate_VASP_structure(structure))

    import tempfile
    import subprocess
    import os
    from subprocess import PIPE

    force_field = 'mm3'
    temp_file_name = tempfile.gettempdir() + '/tinker_temp'+ '_' + str(os.getpid())


    tinker_input_file = open(temp_file_name + '.txyz',mode='w')
    tinker_input_file.write(generate_tinker_txyz_file(supercell_wd))

    tinker_key_file = open(temp_file_name + '.key',mode='w')
    tinker_key_file.write(generate_tinker_key_file(supercell_wd))

    tinker_input_file.close()
    tinker_key_file.close()


    print('filename', tinker_input_file)
    tinker_command = './testgrad ' + tinker_input_file.name + \
                     ' ' + force_field + ' Y N N ' + ' -k ' + tinker_key_file.name

    tinker_process = subprocess.Popen(tinker_command, stdin=PIPE, stderr=PIPE, stdout=PIPE, shell=True)
    (output, err) = tinker_process.communicate()
    tinker_process.wait()

    os.unlink(tinker_input_file.name)
    os.unlink(tinker_key_file.name)

    forces = parse_tinker_forces(output)
    print(forces)
    print(forces.shape)
