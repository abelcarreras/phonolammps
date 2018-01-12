[![PyPI version](https://badge.fury.io/py/phonoLAMMPS.svg)](https://badge.fury.io/py/phonoLAMMPS)

phonoLAMMPS
===========
Calculate the 2nd order force constants using phonopy and LAMMPS.

Main features
-------------
- Command line interface (phonopy like style)
- Python API fully compatible with phonopy
- Use of official LAMMPS python interface
- Simple and easy to use

Requirements
------------
- numpy
- phonopy (https://atztogo.github.io/phonopy/)
- LAMMPS python interface (http://lammps.sandia.gov/doc/Section_python.html)

Optional requirements
---------------------
- matplotlib
- seekpath (https://github.com/giovannipizzi/seekpath)


Installation instructions
--------------------------

1) From source code
```
# python setup.py install --user --prefix=
```

2) From PyPI repository

```
# pip install phonoLAMMPS --user
```

For convenience, you may want to copy (or link) the files inside scripts
folder to a location included in $PATH environment variable

Command line interface
----------------------
phonoLAMMPS has a similar interface to phonopy to allow to easily
calculate the 2nd order force constants and generate the unitcell
from a LAMMPS input file in phonopy/VASP format. All outputs 
are ready to use in phonopy calculations.
Also allows features a quick preview of the phonon 
band structure (requires seekpath). 

```
# phonolammps in.lammps --dim 3 3 3 -c POSCAR_unitcell -p
```

Python API 
----------
Simple python API fully compatible with phonopy.

```
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
```

Contact info
---------------------------------------------------------
Abel Carreras
<br>abelcarreras83@gmail.com

Department of Materials Science and Engineering
<br>Kyoto University