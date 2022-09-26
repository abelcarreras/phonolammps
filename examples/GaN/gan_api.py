# Simple API example to calculate the harmonic force constants

from phonolammps import Phonolammps
import numpy as np

lammps_inp = open('in.lammps').read().split('\n')
print(np.array(lammps_inp))
phlammps = Phonolammps(lammps_inp,
                       supercell_matrix=np.diag([3, 3, 3]),
                       primitive_matrix=np.identity(3))

phlammps.optimize_unitcell(energy_tol=0, force_tol=1e-10)
phlammps.plot_phonon_dispersion_bands()

force_constants = phlammps.get_force_constants()

print(force_constants)