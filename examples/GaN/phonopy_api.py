# Integration with phonopy

from phonolammps import Phonolammps
from phonopy import Phonopy

supercell_matrix = [[3, 0, 0],
                    [0, 3, 0],
                    [0, 0, 3]]

phlammps = Phonolammps('in.lammps',
                       supercell_matrix=supercell_matrix)

unitcell = phlammps.get_unitcell()

force_constants = phlammps.get_force_constants().get_array()

phonon = Phonopy(unitcell,
                 supercell_matrix)

phonon.set_force_constants(force_constants)
phonon.set_mesh([20, 20, 20], is_eigenvectors=True)

phonon.set_total_DOS()
phonon.plot_total_DOS().show()

phonon.set_thermal_properties()
phonon.plot_thermal_properties().show()