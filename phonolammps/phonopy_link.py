from phonopy.api_phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.file_IO import parse_BORN
import numpy as np


# Subclassing PhonopyAtoms to include connectivity & atom types (for tinker & gromacs)
class PhonopyAtomsConnect(PhonopyAtoms):
    def __init__(self, **kwargs):
        # Extract connectivity & atom types
        connectivity = kwargs.pop('connectivity', None)
        atom_types = kwargs.pop('atom_types', None)
        self._connectivity = connectivity
        self._atom_types = atom_types

        super(PhonopyAtomsConnect, self).__init__(**kwargs)

    def get_connectivity(self):
        return self._connectivity

    def get_atom_types(self):
        return self._atom_types


class ForceConstants:
    """
    Define Force constants object
    """
    def __init__(self, force_constants, supercell=np.identity(3)):
        """
        Initialize force constants

        :param force_constants: array matrix containing the force constants (phonopy format)
        :param supercell: 3x3 array (or list of lists) containing the supercell definition
        """

        self._force_constants = np.array(force_constants)
        self._supercell = np.array(supercell)

    def get_array(self):
        """
        get the force constants array in phonopy format

        :return: force constants array
        """
        return self._force_constants

    def get_supercell(self):
        """
        get the supercell (respect to the unit cell) in which the force constants are defined

        :return: 3x3 array containing the supercell
        """
        return self._supercell


def get_phonon(structure,
               NAC=False,
               setup_forces=True,
               super_cell_phonon=np.identity(3),
               primitive_matrix=np.identity(3),
               symmetrize=True):
    """
    Return a phonopy phonon object (instance of the class Phonon)
    
    :param structure: unit cell matrix (lattice vectors in rows)
    :param NAC: (Bool) activate/deactivate Non-analytic corrections
    :param setup_forces: (Bool) decide if pre-calculate harmonic forces in phonon object
    :param super_cell_phonon: 3x3 array containing the supercell to be used to calculate the force constants
    :param primitive_matrix: 3x3 array containing the primitive axis (in rows) which define the primitive cell
    :param symmetrize: decide if symmetrize the force constants
    :return: phonopy phonon object
    """

    phonon = Phonopy(structure, super_cell_phonon,
                     primitive_matrix=primitive_matrix,
                     symprec=1e-5, is_symmetry=symmetrize)

    # Non Analytical Corrections (NAC) from Phonopy [Frequencies only, eigenvectors no affected by this option]

    if setup_forces:
        if structure.get_force_constants() is not None:
            phonon.set_force_constants(structure.get_force_constants().get_array())
        elif structure.get_force_sets() is not None:
            phonon.set_displacement_dataset(structure.get_force_sets().get_dict())
            phonon.produce_force_constants()
            structure.set_force_constants(ForceConstants(phonon.get_force_constants(),
                                                         supercell=structure.get_force_sets().get_supercell()))
        else:
            print('No force sets/constants available!')
            exit()

    if NAC:
        print("Using non-analytical corrections")
        primitive = phonon.get_primitive()
        try:
            nac_params = parse_BORN(primitive, is_symmetry=True)
            phonon.set_nac_params(nac_params=nac_params)
        except OSError:
            print('Required BORN file not found!')
            exit()

    return phonon


def obtain_phonon_dispersion_bands(structure, bands_ranges, force_constants, supercell,
                                   NAC=False, band_resolution=30, band_connection=False,
                                   primitive_matrix=np.identity(3)):
    """
    Get the phonon dispersion bands in phonopy format

    :param structure: unit cell matrix (lattice vectors in rows)
    :param bands_ranges: define the path in the reciprocal space (phonopy format)
    :param force_constants: force constants array ( in phonopy format)
    :param supercell: 3x3 array containing the supercell to be used to calculate the force constants
    :param NAC: (Bool) activate/deactivate Non-analytic corrections
    :param band_resolution: define number of points in path in the reciprocal space
    :param band_connection: decide if bands will be all connected or in segments
    :param primitive_matrix: 3x3 array containing the primitive axis (in rows) which define the primitive cell
    :return:
    """
    phonon = get_phonon(structure, NAC=NAC, setup_forces=False,
                        super_cell_phonon=supercell,
                        primitive_matrix=primitive_matrix)

    phonon.set_force_constants(force_constants)

    bands =[]
    for q_start, q_end in bands_ranges:
        band = []
        for i in range(band_resolution+1):
            band.append(np.array(q_start) + (np.array(q_end) - np.array(q_start)) / band_resolution * i)
        bands.append(band)

    try:
        phonon.run_band_structure(bands, is_band_connection=band_connection, with_eigenvectors=True)
        bands_dict = phonon.get_band_structure_dict()
        bands_phonopy = (bands_dict['qpoints'],
                         bands_dict['distances'],
                         bands_dict['frequencies'],
                         bands_dict['eigenvectors'])

    except AttributeError:
        # phonopy 1.9.x+ support
        phonon.set_band_structure(bands, is_band_connection=band_connection, is_eigenvectors=True)
        bands_phonopy = phonon.get_band_structure()

    return bands_phonopy


def get_primitive_structure(structure, primitive_matrix=np.eye(3)):
    from phonopy.structure.cells import get_primitive
    return get_primitive(structure, primitive_matrix)


def standarize_structure(structure):
    from phonopy.structure.spglib import standardize_cell

    lattice, positions, numbers = standardize_cell(structure, to_primitive=False, no_idealize=False, symprec=1e-5)

    return PhonopyAtoms(positions=positions,
                        numbers=numbers,
                        cell=lattice)
