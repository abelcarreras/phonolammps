.. highlight:: rst

Troubleshooting
===============

Check LAMMPS python API
-----------------------
If there is some problem with LAMMPS installation this is a simple script that
can help you to find it out. Run this script inside one of the example folders
(where in.lammps file is placed) ::

    from lammps import lammps

    lmp1 = lammps()
    lmp1.file("in.lammps")
    lmp1.close()

if everything works as expected you should get an output like this ::

    LAMMPS (15 May 2019)
    OMP_NUM_THREADS environment is not set. Defaulting to 1 thread. (src/comm.cpp:88)
      using 1 OpenMP thread(s) per MPI task
    Lattice spacing in x,y,z = 4.785 2.76262 5.189
    Created triclinic box = (0 0 0) to (3.19 2.76262 5.189) with tilt (-1.595 0 0)
    WARNING: Triclinic box skew is large (src/domain.cpp:194)
      1 by 1 by 1 MPI processor grid
    Created 4 atoms
      create_atoms CPU = 0.00049852 secs
    Reading potential file GaN.tersoff with DATE: 2007-10-22
    Total wall time: 0:00:00


otherwise there may be some trouble with LAMMPS python interface. Check LAMMPS
manual for further information.


Check LAMMPS calculations log
-----------------------------

By default LAMMPS logs are deactivated and not shown during the calculation. If issues appear it may be
useful to check LAMMPS force calculations logs. This is done by using **-logshow** flag. Ex: ::

    $ phonolammps in.lammps --dim 2 2 2 --logshow

