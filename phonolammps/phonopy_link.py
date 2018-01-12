from phonopy.api_phonopy import Phonopy
#from phonopy.structure.atoms import Atoms as PhonopyAtoms
from phonopy.file_IO import parse_BORN
import numpy as np


class ForceConstants:
    def __init__(self, force_constants, supercell=np.identity(3)):

        self._force_constants = np.array(force_constants)
        self._supercell = np.array(supercell)

    def get_array(self):
        return self._force_constants

    def get_supercell(self):
        return self._supercell


def get_phonon(structure,
               NAC=False,
               setup_forces=True,
               super_cell_phonon=np.identity(3),
               primitive_axis=np.identity(3)):

    phonon = Phonopy(structure, super_cell_phonon,
                     primitive_matrix=primitive_axis,
                     symprec=1e-5)

    # Non Analytical Corrections (NAC) from Phonopy [Frequencies only, eigenvectors no affected by this option]

    if setup_forces:
        if structure.get_force_constants() is not None:
            phonon.set_force_constants(structure.get_force_constants().get_array())
        elif structure.get_force_sets() is not None:
            phonon.set_displacement_dataset(structure.get_force_sets().get_dict())
            phonon.produce_force_constants(computation_algorithm="svd")
            structure.set_force_constants(ForceConstants(phonon.get_force_constants(),
                                                         supercell=structure.get_force_sets().get_supercell()))
        else:
            print('No force sets/constants available!')
            exit()

    if NAC:
        print("Using non-analytical corrections")
        primitive = phonon.get_primitive()
        nac_params = parse_BORN(primitive, is_symmetry=True)
        phonon.set_nac_params(nac_params=nac_params)

    return phonon


def obtain_phonon_dispersion_bands(structure, bands_ranges, force_constants, supercell,
                                   NAC=False, band_resolution=30, band_connection=False,
                                   primitive_matrix=np.identity(3)):

    phonon = get_phonon(structure, NAC=NAC, setup_forces=False,
                        super_cell_phonon=supercell,
                        primitive_axis=primitive_matrix)

    phonon.set_force_constants(force_constants)

    bands =[]
    for q_start, q_end in bands_ranges:
        band = []
        for i in range(band_resolution+1):
            band.append(np.array(q_start) + (np.array(q_end) - np.array(q_start)) / band_resolution * i)
        bands.append(band)
    phonon.set_band_structure(bands, is_band_connection=band_connection, is_eigenvectors=True)

    return phonon.get_band_structure()
