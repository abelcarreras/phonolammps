from phonolammps import Phonolammps
from phonolammps.iofile import read_from_file_structure_poscar
import numpy as np

bulk = read_from_file_structure_poscar('POSCAR_unitcell')
phlammps = Phonolammps(bulk,
                       'in.lammps',
                       supercell_matrix=np.diag([3, 3, 3]))
phlammps.plot_phonon_dispersion_bands()
