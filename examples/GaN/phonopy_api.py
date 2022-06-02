# Integration with phonopy

from phonolammps import Phonolammps
from phonopy import Phonopy

phlammps = Phonolammps('in.lammps',
                       supercell_matrix=[[3, 0, 0],
                                         [0, 3, 0],
                                         [0, 0, 3]])

unitcell = phlammps.get_unitcell()
force_constants = phlammps.get_force_constants()
supercell_matrix = phlammps.get_supercell_matrix()

phonon = Phonopy(unitcell,
                 supercell_matrix)

phonon.set_force_constants(force_constants)
phonon.run_mesh([20, 20, 20])

phonon.run_total_dos()
phonon.plot_total_dos()

phonon.run_thermal_properties()
phonon.plot_thermal_properties().show()