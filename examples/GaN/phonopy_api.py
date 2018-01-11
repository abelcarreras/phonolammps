# Integration with phonopy

from phonolammps import Phonolammps
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

cell = [[ 3.1900000572, 0,           0],
        [-1.5950000286, 2.762621076, 0],
        [ 0.0,          0,           5.1890001297]]

scaled_positions=[(0.6666669,  0.3333334,  0.0000000),
                  (0.3333331,  0.6666663,  0.5000000),
                  (0.6666669,  0.3333334,  0.3750000),
                  (0.3333331,  0.6666663,  0.8750000)]

symbols=['Ga', 'Ga', 'N', 'N']

supercell_matrix = [[3, 0, 0],
                    [0, 3, 0],
                    [0, 0, 3]]

unitcell = PhonopyAtoms(symbols=symbols,
                        cell=cell,
                        scaled_positions=scaled_positions)

phlammps = Phonolammps(unitcell,
                       'in.lammps',
                       supercell_matrix=supercell_matrix)

force_constants = phlammps.get_force_constants().get_array()

phonon = Phonopy(unitcell,
                 supercell_matrix)

phonon.set_force_constants(force_constants)
phonon.set_mesh([20, 20, 20], is_eigenvectors=True)

phonon.set_total_DOS()
phonon.plot_total_DOS().show()

phonon.set_thermal_properties()
phonon.plot_thermal_properties().show()