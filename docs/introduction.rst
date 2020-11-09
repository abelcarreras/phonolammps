.. highlight:: rst

Introduction
============

PhonoLAMMPS is a python software designed to interface between *LAMMPS* and *phonopy*. With this software allows
to calculate the 2nd order interatomic force constants with phonopy using the forces obtained by LAMMPS.
For this purpose phonoLAMMPS uses the official LAMMPS python API to link both LAMMPS & phonopy.
This is done in a clean way without generating temporal intermediate files on the disk.

PhonoLAMMPS can be used either as a python module with a very similar interface to phonopy or
via a command line interface using a provided general python script written using argparse.

PhonoLAMMPS has been tested with the following LAMMPS versions:
- 16 March 2018
- 15 May 2019
- 29 Oct 2020.

Main features
-------------
- Command line interface (phonopy like style)
- Python API fully compatible with phonopy
- Use of official LAMMPS python interface
- Simple and easy to use
- New! Finite temperature force constants using DynaPhoPy

