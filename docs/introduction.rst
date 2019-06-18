.. highlight:: rst

Introduction
============

phonoLAMMPS is a python software designed to interface between LAMMPS and phonopy to calculate
2n order interatomic force constants using phonopy with the forces calculated by LAMMPS. For this purpose phonoLAMMPS
uses the officially supported LAMMPS python API to link both LAMMPS & phonopy. This is done in a clean way without
intermediate files.

PhonoLAMMPS can be used either as a python module with a very similar inteface to phonopy or
via a command line interface using a python script.

Main features
-------------
- Command line interface (phonopy like style)
- Python API fully compatible with phonopy
- Use of official LAMMPS python interface
- Simple and easy to use
- New! Finite temperature force constants using DynaPhoPy

