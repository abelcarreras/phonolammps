.. highlight:: rst

Get started
===========

Python API
----------
phonoLAMMPS has been designed to be very similar to phonopy python interface.
The different objects generated are fully compatible with phonopy.
The procedure to obtain the force constants from LAMMPS forces is the following:

Loading phonolammps module ::

    from phonolammps import Phonolammps

Creating instance of main phonoLAMMPS class. I this call you have to introduce a lammps
input that contains the definition of the crystal unit cell the unitcell and potentials.
For this you have two options:
1) Use a LAMMPS input file (in.lammps) ::

    phlammps = Phonolammps('in.lammps',
                           supercell_matrix=[[2, 0, 0],
                                             [0, 2, 0],
                                             [0, 0, 2]],
                           primitive_matrix=[[0.0, 0.5 ,0.5],
                                             [0.5, 0.0, 0.5],
                                             [0.5, 0.5, 0.0])

2) Use a python list containing the lammps commands (one command per line) ::

    list_of_commands = open('in.lammps').read().split('\n')
    phlammps = Phonolammps(list_of_commands,
                           supercell_matrix=[[2, 0, 0],
                                             [0, 2, 0],
                                             [0, 0, 2]])
                           primitive_matrix=[[0.0, 0.5 ,0.5],
                                             [0.5, 0.0, 0.5],
                                             [0.5, 0.5, 0.0])

This second option is convenient if you want to edit lammps commands in python scripting.
Also in the creation of the instance of this class you have to define the supercell expansion
to be used in the force constants calculation (ex: 2x2x2), and if necessary, the primitive cell
axis matrix. If not defined by default will be the identity matrix (primitive cell = unit cell).


Get the data needed for phonopy ::

    unitcell = phlammps.get_unitcell()
    force_constants = phlammps.get_force_constants()
    supercell_matrix = phlammps.get_supercell_matrix()


From this you have all the information you need for phonopy calculations ::

    from phonopy import Phonopy
    phonon = Phonopy(unitcell,
                     supercell_matrix)

    phonon.set_force_constants(force_constants)
    phonon.set_mesh([20, 20, 20])

    phonon.set_total_DOS()
    phonon.plot_total_DOS().show()

    phonon.set_thermal_properties()
    phonon.plot_thermal_properties().show()


Command line interface
----------------------
To use phonoLAMMPS from command line interface you need a LAMMPS input file containing the
definition of the unit cell and potentials.

example: ::

    units           metal

    boundary        p p p

    box tilt large

    atom_style      atomic

    read_data       data.si

    pair_style      tersoff
    pair_coeff      * * SiCGe.tersoff  Si(C)

    neighbor	0.3 bin

*Notice that run command should not be written in the input file. **read data** command can be used to define
the atoms coordinates in another file.


Phonolammps script uses argparse to provide a clean command line interface using flags
The different options available are displayed by using **-h** flag ::

    $ phonolammps -h
    usage: phonolammps [-h] [-o file] [--dim N N N] [-p] [-c file]
                       [-pa F F F F F F F F F] [-t F] [--amplitude F]
                       [--total_time F] [--relaxation_time F] [--timestep F]
                       [--logshow] [--no_symmetrize] [--use_NAC]
                       [--write_force_sets]
                       lammps_file

    phonoLAMMPS options

    positional arguments:
      lammps_file           lammps input file

    optional arguments:
      -h, --help            show this help message and exit
      -o file               force constants output file [default: FORCE_CONSTANTS]
      --dim N N N           dimensions of the supercell
      -p                    plot phonon band structure
      -c file, --cell file  generates a POSCAR type file containing the unit cell
      -pa F F F F F F F F F, --primitive_axis F F F F F F F F F
                            primitive axis
      -t F                  temperature in K
      --amplitude F         displacement distance [default: 0.01 angstrom]
      --total_time F        total MD time in picoseconds [default: 20 ps]
      --relaxation_time F   MD relaxation time in picoseconds [default: 5 ps]
      --timestep F          MD time step in picoseconds [default: 0.001 ps]
      --logshow             show LAMMPS & dynaphopy log on screen
      --no_symmetrize       deactivate force constant symmetrization
      --use_NAC             include non analytical corrections (Requires BORN file
                            in work directory)
      --write_force_sets    write FORCE_SETS file


A simple example for crystalline silicon using a 2x2x2 supercell would be ::

    phonolammps in.lammps --dim 2 2 2 -pa 0.0 0.5 0.5 0.5 0.0 0.5 0.5 0.5 0.0 -c POSCAR_unitcell -p

where **in.lammps** is a lammps input containing the unit cell, --pa are the primitive axis in matrix
format written in one line and POSCAR_unitcell is the unitcell (the same written in lammps input) written
in VASP format to be used in phonopy calculations.
The output of this script is a file named **FORCE_CONSTANTS** that contains the interatomic 2n order
force constants in phonopy format.

