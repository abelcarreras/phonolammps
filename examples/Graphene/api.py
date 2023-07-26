# API example to calculate finite temperature (300K) force constants
# using dynaphopy (quasiparticle theory)

from phonolammps import Phonolammps
import numpy as np


sm = np.diag([10, 10, 1])
pm = np.identity(3)

# calculate harmonic force constants with phonoLAMMPS
phlammps = Phonolammps('in.lammps',
                       supercell_matrix=sm,
                       primitive_matrix=pm)


bl = {
    "ranges": [[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
               [[0.5, 0.0, 0.0], [0.3333333333, 0.3333333333333, 0.0]],
               [[0.3333333333, 0.3333333333333, 0.0], [0.0, 0.0, 0.0]]],
    "labels": [("GAMMA", "M"), ("M", "K"), ("K", "GAMMA")]
}
phlammps.plot_phonon_dispersion_bands(bands_and_labels=bl)
