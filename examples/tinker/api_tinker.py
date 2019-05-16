# API example to calculate force constants with tinker (on development)
# Only works for Orthorhombic crystals

from phonolammps import PhonoTinker
import numpy as np

phtinker = PhonoTinker(txyz_input_file='structure.txyz',
                       key_input_file='structure.key',
                       force_field_file='mm3.prm',
                       supercell_matrix=np.diag([2, 2, 2]),
                       show_progress=True)

phtinker.write_force_constants()
print('Force constants generated')

phtinker.write_unitcell_POSCAR()
print('POSCAR generated')

#phtinker.plot_phonon_dispersion_bands()

