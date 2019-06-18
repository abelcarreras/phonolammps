import os
import numpy as np

from phonopy.structure.atoms import PhonopyAtoms
from phonolammps.phonopy_link import PhonopyAtomsTinker


def mass_to_symbol(mass, tolerance=5e-1):
    from phonopy.structure.atoms import atom_data

    for element in atom_data:
        if element[3] is not None and abs(mass - element[3]) < tolerance:
            return element[1]

    return 'H'  # in case no match found use H as wildcard


def get_structure_from_poscar(file_name, number_of_dimensions=3):
    """
    Read crystal structure from a VASP POSCAR type file

    :param file_name: POSCAR filename
    :param number_of_dimensions: number of dimensions of the crystal structure
    :return: Atoms (phonopy) type object containing the crystal structure
    """
    # Check file exists
    if not os.path.isfile(file_name):
        print('Structure file does not exist!')
        exit()

    # Read from VASP POSCAR file
    poscar_file = open(file_name, 'r')
    data_lines = poscar_file.read().split('\n')
    poscar_file.close()

    multiply = float(data_lines[1])
    direct_cell = np.array([data_lines[i].split()
                            for i in range(2, 2+number_of_dimensions)], dtype=float)
    direct_cell *= multiply
    scaled_positions = None
    positions = None

    try:
        number_of_types = np.array(data_lines[3+number_of_dimensions].split(),dtype=int)

        coordinates_type = data_lines[4+number_of_dimensions][0]
        if coordinates_type == 'D' or coordinates_type == 'd' :

            scaled_positions = np.array([data_lines[8+k].split()[0:3]
                                         for k in range(np.sum(number_of_types))],dtype=float)
        else:
            positions = np.array([data_lines[8+k].split()[0:3]
                                  for k in range(np.sum(number_of_types))],dtype=float)

        atomic_types = []
        for i,j in enumerate(data_lines[5].split()):
            atomic_types.append([j]*number_of_types[i])
        atomic_types = [item for sublist in atomic_types for item in sublist]

    # Old style POSCAR format
    except ValueError:
        number_of_types = np.array(data_lines[5].split(), dtype=int)
        coordinates_type = data_lines[6][0]
        if coordinates_type == 'D' or coordinates_type == 'd':
            scaled_positions = np.array([data_lines[7+k].split()[0:3]
                                         for k in range(np.sum(number_of_types))], dtype=float)
        else:
            positions = np.array([data_lines[7+k].split()[0:3]
                                  for k in range(np.sum(number_of_types))], dtype=float)

        atomic_types = []
        for i,j in enumerate(data_lines[0].split()):
            atomic_types.append([j]*number_of_types[i])
        atomic_types = [item for sublist in atomic_types for item in sublist]

    return PhonopyAtoms(symbols=atomic_types,
                        scaled_positions=scaled_positions,
                        cell=direct_cell)


def get_structure_from_lammps(command_list, show_log=False):
    """
    Get the crystal structure from lammps input

    :param command_list: LAMMPS input commands in list (one item for line)
    :return: numpy array matrix with forces of atoms [Natoms x 3]
    """
    from lammps import lammps

    cmd_list = ['-log', 'none']
    if not show_log:
        cmd_list += ['-echo', 'none', '-screen', 'none']

    lmp = lammps(cmdargs=cmd_list)

    lmp.commands_list(command_list)
    lmp.command('run 0')

    na = lmp.get_natoms()

    xlo =lmp.extract_global("boxxlo", 1)
    xhi =lmp.extract_global("boxxhi", 1)
    ylo =lmp.extract_global("boxylo", 1)
    yhi =lmp.extract_global("boxyhi", 1)
    zlo =lmp.extract_global("boxzlo", 1)
    zhi =lmp.extract_global("boxzhi", 1)
    xy =lmp.extract_global("xy", 1)
    yz =lmp.extract_global("yz", 1)
    xz =lmp.extract_global("xz", 1)

    unitcell = np.array([[xhi-xlo, xy,  xz],
                           [0,  yhi-ylo,  yz],
                           [0,   0,  zhi-zlo]]).T

    positions = lmp.gather_atoms("x", 1, 3)
#    type_mass = lmp.gather_atoms("mass", 1, 1)
    type_mass = lmp.extract_atom("mass", 2)

    type = lmp.gather_atoms("type", 0, 1)

    positions = np.array([positions[i] for i in range(na * 3)]).reshape((na, 3))
    masses = np.array([type_mass[type[i]] for i in range(na)])
    symbols = [mass_to_symbol(masses[i]) for i in range(na)]

    return PhonopyAtoms(positions=positions,
                        masses=masses,
                        symbols=symbols,
                        cell=unitcell)

def generate_VASP_structure(structure, scaled=True, supercell=(1, 1, 1)):

    cell = structure.get_cell()
    types = structure.get_chemical_symbols()

    elements = [types[0]]
    elements_count = [1]

    for t in types[1:]:
        if t == elements[-1]:
            elements_count[-1] += 1
        else:
            elements.append(t)
            elements_count.append(1)

    #atom_type_unique = np.unique(types, return_index=True)

    # To use unique without sorting
    #sort_index = np.argsort(atom_type_unique[1])
    #elements = np.array(atom_type_unique[0])[sort_index]

    #elements_count= np.diff(np.append(np.array(atom_type_unique[1])[sort_index], [len(types)]))
    #print(elements_count)

    vasp_POSCAR = 'Generated using phonoLAMMPS\n'
    vasp_POSCAR += '1.0\n'
    for row in cell:
        vasp_POSCAR += '{0:20.10f} {1:20.10f} {2:20.10f}\n'.format(*row)
    vasp_POSCAR += ' '.join(elements)
    vasp_POSCAR += ' \n'
    vasp_POSCAR += ' '.join([str(i) for i in elements_count])

    if scaled:
        scaled_positions = structure.get_scaled_positions()
        vasp_POSCAR += '\nDirect\n'
        for row in scaled_positions:
            vasp_POSCAR += '{0:15.15f}   {1:15.15f}   {2:15.15f}\n'.format(*row)

    else:
        positions = structure.get_positions()
        vasp_POSCAR += '\nCartesian\n'
        for row in positions:
            vasp_POSCAR += '{0:20.10f} {1:20.10f} {2:20.10f}\n'.format(*row)

    return vasp_POSCAR


def get_structure_from_txyz(file_name, key_file):

    tinker_file = open(file_name, 'r')

    coordinates = []
    atomic_numbers = []
    atomic_elements = []
    atom_types = []
    connectivity = []

    number_of_atoms = int(tinker_file.readline().split()[0])

    # print(number_of_atoms)
    for i in range(number_of_atoms):
        line = tinker_file.readline().split()
        # Check if txyz contains cell parameters
        if line[1].replace('.','',1).isdigit():
            line = tinker_file.readline().split()

        coordinates.append(line[2:5])
        atomic_numbers.append(int(line[0]))
        atomic_elements.append(line[1])
        atom_types.append(line[5])
        connectivity.append([int(f) for f in line[6:]])
    # print(np.array(coordinates,dtype=float))

    tinker_file.close()

    lines = open(key_file, 'r').readlines()
    # print lines

    params = []
    for label in ['a-axis', 'b-axis', 'c-axis', 'alpha', 'beta', 'gamma']:
        for line in lines:
            if label in line.lower():
                params.append(float(line.split()[1]))
    a, b, c, alpha, beta, gamma = params

    a1 = a
    b1 = b*np.cos(np.deg2rad(gamma))
    b2 = np.sqrt(b**2 - b1**2)
    c1 = c*np.cos(np.deg2rad(beta))
    c2 = (b*c*np.cos(np.deg2rad(alpha)) - b1*c1)/b2
    c3 = np.sqrt(c**2-c1**2-c2**2)

    unitcell = [[a1, 0,  0],
                [b1, b2, 0],
                [c1, c2, c3]]

    return PhonopyAtomsTinker(positions=np.array(coordinates, dtype=float),
                              symbols=atomic_elements,
                              cell=unitcell,
                              connectivity=connectivity,
                              atom_types=atom_types)


def generate_tinker_txyz_file(structure):
    pos = structure.get_positions()
    sym = structure.get_chemical_symbols()
    con = structure.get_connectivity()
    typ = structure.get_atom_types()

    tinker_txt = '{}      MM2 parameters\n'.format(len(pos))
    for i, p in enumerate(pos):
        tinker_txt += ' {} {} '.format(i+1, sym[i]) + '{} {} {} '.format(*p) + ' {} '.format(typ[i]) + ' '.join(np.array(con[i], dtype=str)) + '\n'

    return tinker_txt


def generate_tinker_key_file(structure,
                             archive='overwrite',
                             wdv_cutoff=20):

    cell = structure.get_cell()

    a, b, c = np.sqrt(np.dot(cell, cell.transpose()).diagonal())
    alpha = np.rad2deg(np.arccos(np.vdot(cell[1], cell[2]) / b / c))
    beta = np.rad2deg(np.arccos(np.vdot(cell[2], cell[0]) / c / a))
    gamma = np.rad2deg(np.arccos(np.vdot(cell[0], cell[1]) / a / b))

    tinker_txt = '# Cell parameters\n'
    tinker_txt += 'a-axis  {}\n'.format(a)
    tinker_txt += 'b-axis  {}\n'.format(b)
    tinker_txt += 'c-axis  {}\n'.format(c)

    tinker_txt += 'ALPHA  {}\n'.format(alpha)
    tinker_txt += 'BETA   {}\n'.format(beta)
    tinker_txt += 'GAMMA  {}\n'.format(gamma)

    tinker_txt += '# Other parameters\n'.format(c)
    tinker_txt += 'VDW-CUTOFF {}\n'.format(wdv_cutoff)
    tinker_txt += 'archive {}\n'.format(archive)

    return tinker_txt


def parse_tinker_forces(output):

    data = output.decode().split('Anlyt')[1:-2]

    gradient = []
    for d in data:
        gradient.append(d.split()[1:4])
    forces = -np.array(gradient, dtype=float)
    return forces
