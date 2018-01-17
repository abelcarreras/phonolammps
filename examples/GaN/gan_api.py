# Simple API example

from phonolammps import Phonolammps
import numpy as np

phlammps = Phonolammps('in.lammps',
                       supercell_matrix=np.diag([3, 3, 3]),
                       primitive_matrix=np.identity(3))

phlammps.plot_phonon_dispersion_bands()

force_constants = phlammps.get_force_constants()

print force_constants