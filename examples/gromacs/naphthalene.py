from phonolammps import PhonoGromacs
import numpy as np


gmx_params = {'rvdw': 0.28,      # nm
              'rlist': 0.28,     # nm
              'rcoulomb': 0.28}  # nm

phg = PhonoGromacs('unitcell_whole.g96',
                   supercell_matrix=np.identity(3) * 4,
                   displacement_distance=0.15, # angs
                   gmx_params=gmx_params,
                   itp_file='params.itp',
                   omp_num_threads=6,
                   show_progress=True)

phg.write_unitcell_POSCAR()
phg.plot_phonon_dispersion_bands()
phg.write_force_constants()
phonon = phg.get_phonopy_phonon()
phonon.run_mesh([40, 40, 40])
phonon.run_total_dos()
phonon.plot_total_dos().show()
