import numpy as np


def diff_matrix(array_1, array_2, cell_size):
    """
    :param array_1: supercell scaled positions respect unit cell
    :param array_2: supercell scaled positions respect unit cell
    :param cell_size: diference between arrays accounting for periodicity
    :return:
    """
    array_1_norm = np.array(array_1) / np.array(cell_size, dtype=float)[None,:]
    array_2_norm = np.array(array_2) / np.array(cell_size, dtype=float)[None,:]

    return array_2_norm - array_1_norm


def phonopy_order(i, size):
    x = np.mod(i, size[0])
    y = np.mod(i, size[0]*size[1])//size[1]
    z = np.mod(i, size[0]*size[1]*size[2])//(size[1]*size[0])
    k = i//(size[1]*size[0]*size[2])

    return np.array([x, y, z, k])


def get_correct_arrangement(reference, structure, supercell_matrix):

    # print structure.get_scaled_positions()
    scaled_coordinates = []
    for coordinate in reference:
        trans = np.dot(coordinate, np.linalg.inv(structure.get_cell()))
        scaled_coordinates.append(np.array(trans.real, dtype=float))

    number_of_cell_atoms = structure.get_number_of_atoms()
    number_of_supercell_atoms = len(scaled_coordinates)
    supercell_dim = np.diag(supercell_matrix)

    # print 'atom', number_of_cell_atoms, number_of_supercell_atoms
    unit_cell_scaled_coordinates = scaled_coordinates - np.array(scaled_coordinates, dtype=int)

    # Map supercell atoms to unitcell
    atom_unit_cell_index = []
    for coordinate in unit_cell_scaled_coordinates:
        # Only works for non symmetric cell (must be changed)

        diff = np.abs(np.array([coordinate]*number_of_cell_atoms) - structure.get_scaled_positions())

        diff[diff >= 0.5] -= 1.0
        diff[diff < -0.5] += 1.0

        index = np.argmin(np.linalg.norm(diff, axis=1))

        atom_unit_cell_index.append(index)
    atom_unit_cell_index = np.array(atom_unit_cell_index)

    original_conf = np.array([phonopy_order(j, supercell_dim)[:3] for j in range(number_of_supercell_atoms)])

    template = []
    lp_coordinates = []
    for i, coordinate in enumerate(scaled_coordinates):
        lattice_points_coordinates = coordinate - structure.get_scaled_positions()[atom_unit_cell_index[i]]

        for k in range(3):
            if lattice_points_coordinates[k] > supercell_dim[k] - 0.5:
                lattice_points_coordinates[k] = lattice_points_coordinates[k] - supercell_dim[k]
            if lattice_points_coordinates[k] < -0.5:
                lattice_points_coordinates[k] = lattice_points_coordinates[k] + supercell_dim[k]

        comparison_cell = np.array([lattice_points_coordinates]*number_of_supercell_atoms)
        diference = np.linalg.norm(diff_matrix(original_conf, comparison_cell, supercell_dim), axis=1)
        template.append(np.argmin(diference) + atom_unit_cell_index[i]*number_of_supercell_atoms//number_of_cell_atoms)

        lp_coordinates.append(lattice_points_coordinates)
    template = np.array(template)

    if len(np.unique(template)) < len(template):
        print ('Something wrong with crystal structure!\n'
               'POSCAR & LAMMPS structure do not match')
        print ('unique: {} / {}'.format(len(np.unique(template)), len(template)))
        exit()

    return template
