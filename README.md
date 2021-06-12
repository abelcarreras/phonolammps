[![PyPI version](https://badge.fury.io/py/phonoLAMMPS.svg)](https://badge.fury.io/py/phonoLAMMPS)
[![Downloads](http://pepy.tech/badge/phonolammps)](http://pepy.tech/project/phonolammps)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3940625.svg)](https://doi.org/10.5281/zenodo.3940625)


phonoLAMMPS
===========
Calculate the harmonic interatomic force constants using phonopy and LAMMPS.  
Online manual: https://phonolammps.readthedocs.io

Main features
-------------
- Command line interface (phonopy like style)
- Python API fully compatible with phonopy
- Use of official LAMMPS python interface
- Simple and easy to use
- Finite temperature force constants using DynaPhoPy

Requirements
------------
- python 2.7.x/3.x
- numpy
- phonopy (https://atztogo.github.io/phonopy/)
- LAMMPS python interface (https://lammps.sandia.gov/doc/Python_library.html)

Optional requirements for phonon band structure preview
-------------------------------------------------------
- matplotlib
- seekpath (https://github.com/giovannipizzi/seekpath)

Optional requirements for finite temperature FC calculations
------------------------------------------------------------
- DynaPhoPy (https://github.com/abelcarreras/DynaPhoPy)

Optional requirements for tinker
------------------------------------------------------------
- Tinker testgrad (https://dasher.wustl.edu/tinker/)


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
calculate the 2nd order force constants and generate the crystal unitcell
from a LAMMPS input file in VASP/POSCAR format. All outputs
are fully compatible and ready to use in phonopy calculations.
Also features a quick preview of the phonon
band structure (requires seekpath). 

```
# phonolammps in.lammps --dim 2 2 2 -c POSCAR_unitcell -p
```
Additionally phonoLAMMPS allows to easily calculate finite temperature force constants 
from molecular dynamics by quasiparticle theory (requires dynaphopy).
```
# phonolammps in.lammps --dim 2 2 2  -c POSCAR_unitcell -p -t 300       (at 300 K)
```
The obtained *FORCE_CONSTANTS* and *POSCAR_unitcell* can be used in phonopy using 
**--readfc** option for more advanced calculations.
 ```
# phonopy --dim="2 2 2" --readfc -c POSCAR_unitcell band.conf
```

Python API 
----------
Simple python API fully compatible with phonopy.

```
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
phonon.set_mesh([20, 20, 20])

phonon.set_total_DOS()
phonon.plot_total_DOS().show()

phonon.set_thermal_properties()
phonon.plot_thermal_properties().show()
```

Contact info
---------------------------------------------------------
Abel Carreras
<br>abelcarreras83@gmail.com

Donostia International Physics Center (DIPC)
<br>Donostia-San Sebastian (Spain)