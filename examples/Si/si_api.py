# API example to calculate finite temperature (300K) force constants
# using dynaphopy (quasiparticle theory)

from phonolammps import Phonolammps
from dynaphopy.interface.lammps_link import generate_lammps_trajectory
from dynaphopy.interface.phonopy_link import ForceConstants
from dynaphopy import Quasiparticle
from dynaphopy.atoms import Structure

import numpy as np

# calculate harmonic force constants with phonoLAMMPS
phlammps = Phonolammps('in.lammps',
                       supercell_matrix=np.diag([2, 2, 2]),
                       primitive_matrix=[[0.0, 0.5, 0.5],
                                         [0.5, 0.0, 0.5],
                                         [0.5, 0.5, 0.0]],
                       )

phlammps.plot_phonon_dispersion_bands()

force_constants = phlammps.get_force_constants()

# Print harmonic force constants
print('harmonic force constants')
print(force_constants)

structure = phlammps.get_unitcell()

# define structure for dynaphopy
dp_structure = Structure(cell=structure.get_cell(),  # cell_matrix, lattice vectors in rows
                         scaled_positions=structure.get_scaled_positions(),
                         atomic_elements=structure.get_chemical_symbols(),
                         primitive_matrix=phlammps.get_primitve_matrix(),
                         force_constants=ForceConstants(force_constants,
                                                        supercell=phlammps.get_supercell_matrix()))


# calculate trajectory for dynaphopy with lammps
trajectory = generate_lammps_trajectory(dp_structure, 'in.lammps',
                                        total_time=20,      # ps
                                        time_step=0.001,    # ps
                                        relaxation_time=5,  # ps
                                        silent=False,
                                        supercell=[2, 2, 2],
                                        memmap=False,
                                        velocity_only=True,
                                        temperature=300)

# set dynaphopy calculation
calculation = Quasiparticle(trajectory)
calculation.select_power_spectra_algorithm(2)  # select FFT algorithm

calculation.get_renormalized_phonon_dispersion_bands()
renormalized_force_constants = calculation.get_renormalized_force_constants()

# Print phonon band structure
calculation.plot_renormalized_phonon_dispersion_bands()

# Print renormalized force constants
print('renormalized force constants at 300K')
print(renormalized_force_constants)