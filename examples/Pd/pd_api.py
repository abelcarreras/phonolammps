# Example to calculate renormalized force constants at finite temperature
# by using molecular dynamics & phonon quasiparticle approach
# Requires: lammps, phonopy & dynaphopy

from phonolammps import Phonolammps
from dynaphopy import Quasiparticle
from dynaphopy.atoms import Structure
from dynaphopy.interface.lammps_link import generate_lammps_trajectory
from dynaphopy.interface.phonopy_link import ForceConstants
from contextlib import contextmanager
from phonopy.file_IO import write_FORCE_CONSTANTS, write_force_constants_to_hdf5

import numpy as np
import os
import sys


@contextmanager
def silence_stdout():

    new_target = open(os.devnull, "w")
    old_target, sys.stdout = sys.stdout, new_target
    try:
        yield new_target
    finally:
        sys.stdout = old_target

# input parameters
supercell = [2, 2, 2]
primitive_mat = [[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 1]]
temperature = 800

# calculate harmonic force consntants with phonolammps
phlammps = Phonolammps('in.lammps',
                       supercell_matrix=np.diag(supercell),
                       primitive_matrix=primitive_mat)

force_constants = phlammps.get_force_constants()
unitcell = phlammps.get_unitcell()

structure = Structure(cell=unitcell.get_cell(),  # cell_matrix, lattice vectors in rows
                      scaled_positions=unitcell.get_scaled_positions(),
                      atomic_elements=unitcell.get_chemical_symbols(),
                      primitive_matrix=primitive_mat)

structure.set_force_constants(ForceConstants(force_constants,
                                             supercell=np.diag(supercell)))

# generate LAMMPS MD trajectory (requires dynaphopy development)
trajectory = generate_lammps_trajectory(structure, 'in.lammps',
                                        total_time=20,
                                        time_step=0.001,
                                        relaxation_time=5,
                                        silent=False,
                                        supercell=supercell,
                                        memmap=False,
                                        velocity_only=False,
                                        temperature=temperature)

# Calculate renormalized force constants with dynaphopy
calculation = Quasiparticle(trajectory)
with silence_stdout():
    calculation.select_power_spectra_algorithm(2)  # 1: Max. Entrop.  2:FFT
    calculation.plot_renormalized_phonon_dispersion_bands()
    renormalized_force_constants = calculation.get_renormalized_force_constants()

# Save renormalized force constants and unitcell to disk
write_FORCE_CONSTANTS(renormalized_force_constants.get_array(),
                      filename='FORCE_CONSTANTS_{}K'.format(temperature))
phlammps.write_unitcell_POSCAR(filename='POSCAR_unitcell')
