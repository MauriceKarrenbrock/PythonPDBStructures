# -*- coding: utf-8 -*-
# pylint: disable=too-many-branches
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Functions for geometrical treatment of PDB files

There are functions to for example calculate the center of
mass COM of a protein structure
"""

import Bio.PDB
import numpy as np

import PythonPDBStructures.important_lists


def get_center_of_mass(structure, geometric=False):
    """calculate center of mass of a biopython structure

    Given a biopython structure (protein, chain)
    or an iterable of atoms (object with element and
    coord attributes)
    it returns the center of mass COM.
    If geometric = True it calculates the geomatrical center
    (all masses are equal)

    Parameters
    -----------
    biopython_structure : biopython entity or iterable of atoms
        can be a protein, model, chain, residue, iterable of atoms
        an atom is an object with element and coord attributes
        coord must be an iterable with (x,y,z)
    geometric : bool, optional
        False (default) calculate the true center of mass
        True calculates the geometric center (all masses are
        equal)

    Returns
    ----------
    numpy.array
        [x, y, z] x, y, z are float
    """

    atom_weights = PythonPDBStructures.important_lists.atom_weights

    # Structure, Model, Chain, Residue
    if isinstance(structure, Bio.PDB.Entity.Entity):
        atom_list = structure.get_atoms()
    # iterable of Atoms
    elif hasattr(structure, '__iter__'):
        atom_list = structure
    else:  # Some other weirdo object
        raise ValueError(
            'Center of Mass can only be calculated from the following objects:\n'
            'Structure, Model, Chain, Residue, iterable of Atoms.')

    positions = [[], [],
                 []]  # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]
    masses = []

    for atom in atom_list:

        #some times atom element fails with strange atom names
        #So I create a "backup name" by removing numbers from atom.name
        #capitalizing it and taking the first 2 chars
        backup_name = ''.join(i for i in atom.name if not i.isdigit())
        backup_name = backup_name.capitalize()
        if len(backup_name) > 2:
            backup_name = backup_name[0:2]

        if atom.element.capitalize() in atom_weights.keys():
            atom.mass = atom_weights[atom.element.capitalize()]

        elif backup_name in atom_weights.keys():
            atom.mass = atom_weights[backup_name]

        else:
            atom.mass = 'ukn'

        masses.append(atom.mass)

        for i in range(len(atom.coord)):
            positions[i].append(atom.coord[i])

    # If there is a single atom with undefined mass complain loudly.
    if 'ukn' in set(masses) and not geometric:
        raise ValueError("Some Atoms don't have an element assigned.\n" +
                         'Try adding them manually or calculate the ' +
                         'geometrical center of mass instead.')

    if geometric:
        com = [sum(coord_list) / len(masses) for coord_list in positions]

    else:
        w_pos = [[], [], []]
        for atom_index, atom_mass in enumerate(masses):
            w_pos[0].append(positions[0][atom_index] * atom_mass)
            w_pos[1].append(positions[1][atom_index] * atom_mass)
            w_pos[2].append(positions[2][atom_index] * atom_mass)
        com = [sum(coord_list) / sum(masses) for coord_list in w_pos]

    #make it a numpy array
    com = np.array(com)

    return com
