# -*- coding: utf-8 -*-
# pylint: disable=too-many-branches
# pylint: disable=too-many-locals
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
import mdtraj
import numpy as np
from simtk.openmm import unit

import PythonPDBStructures.important_lists as important_lists


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

    atom_weights = important_lists.atom_weights

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

    unknown = []
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

            unknown.append(atom)

        masses.append(atom.mass)

        for i in range(len(atom.coord)):
            positions[i].append(atom.coord[i])

    # If there is a single atom with undefined mass complain loudly.
    if 'ukn' in set(masses) and not geometric:
        raise ValueError("Some Atoms don't have an element assigned.\n" +
                         'Try adding them manually or calculate the ' +
                         f'geometrical center of mass instead.\n{unknown}')

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


def get_nearest_neighbors_residues_with_mdtraj(mdtraj_trajectory,
                                               *,
                                               target_resnum=None,
                                               target_resname=None,
                                               cutoff=4.5 * unit.angstrom,
                                               custom_selection=None):
    """get the nearest reighboring residues

    starting from a mdtraj trajectory it gives you the
    residues that have at least one atom that had a distance lower or equal to
    `cutoff_angstom` from at least one atom of the target residue

    it will only consider the first frame in the trajectory!

    this function is very handy if working with openmm (see examples)

    Parameters
    -------------
    mdtraj_trajectory : mdtraj trajectory
    target_resnum : int, optional
        the residue number of the residue (resSeq) (numbering starts from 1)
        from which you want the nearest
        neighbors, it is suggested to give `target_resnum` as input
        instead of `target_resname`, if `target_resnum` is given `target_resname`
        is ignored
    target_resname : str, optional
        the residue name of the residue from which you want the nearest
        neighbors, it is suggested to give `target_resnum` as input
        instead of `target_resname`, if `target_resnum` is given `target_resname`
        is ignored
    cutoff : openmm.unit.Quantity or float, optional, default=4.5 * unit.angstrom
        the cutoff to decide if an atom is a neighbor or not, if the input is not a
        openmm.unit.Quantity it will be considered in angstrom
    custom_selection : str, optional, default=None
        for default the nearest neighbours will be searched only in the 'protein'
        mdtraj selection, if you give a mdtraj selection string the nearest neighbors
        will be searched in that custom selection

    Returns
    ----------
    list(int)
        the list of residue numbers (resSeq) (1 indexed)

    Raises
    ------------
    ValueError
        if both `target_resname` and `target_resnum` are not given as input

    Notes
    -----------
    If you give `target_resname` and there are more than one residue with this name
    results are unpredictable

    it will only consider the first frame in the trajectory!

    Examples
    --------------
    ```
    import mdtraj

    traj = mdtraj.load('../system_with_h.pdb')

    nn = get_nearest_neighbors_residues_with_mdtraj(traj, target_resname='LIG')

    nn = [0, 5, ...]
    ```

    ```
    top = mdtraj.Topology.from_openmm(openmm_topology)
    traj = mdtraj.Trajectory([openmm_positions / unit.nanometers], top)
    ```
    """

    top = mdtraj_trajectory.topology

    if unit.is_quantity(cutoff):

        cutoff = cutoff.value_in_unit(unit.nanometer)

    else:

        # from angstrom to nanometers
        cutoff /= 10.0

    if target_resnum is not None:

        #from 1 indexed resseq to 0 indexed index
        target_resnum -= 1

    elif target_resname is not None:

        target_resnum = top.select(f'resname {target_resname}')
        target_resnum = top.atom(target_resnum[0]).residue.index

    else:

        raise ValueError('You must  give a target_resname or a target_resnum')

    if custom_selection is None:

        custom_selection = 'protein'

    protein_atoms = top.select(custom_selection)

    protein_res = list(set(top.atom(i).residue.index for i in protein_atoms))

    pairs = [(target_resnum, i) for i in protein_res]

    contacts = mdtraj.compute_contacts(mdtraj_trajectory,
                                       contacts=pairs,
                                       scheme='closest')

    #filter the interesting contacts
    interesting_residues = []
    for dist, residues in zip(contacts[0][0], contacts[1]):

        if dist <= cutoff:

            for i in residues:

                #I am more interested in the pdb 1 based index
                #than in the mdtraj 0 based index
                interesting_residues.append(i + 1)

    #remove possible duplicates
    #and remove possible contacts between the target and its self
    interesting_residues = set(interesting_residues)
    interesting_residues.discard(target_resnum + 1)
    interesting_residues = list(interesting_residues)

    return interesting_residues


def get_nearest_neighbors_residues_with_biopython(structure,
                                                  *,
                                                  target_resnum=None,
                                                  target_resname=None,
                                                  ignore_resnums=None,
                                                  ignore_resnames=None,
                                                  cutoff_angstom=4.5,
                                                  ignore_hetatms=False):
    """get the nearest reighboring residues

    starting from a biopython structure it gives you the
    residues that have at least one atom that had a distance lower or equal to
    `cutoff_angstom` from at least one atom of the target residue

    Parameters
    -------------
    structure : biopython structure
        the structure parsed from a pdb or mmcif file via biopython
    target_resnum : int, optional
        the residue number of the residue from which you want the nearest
        neighbors, it is suggested to give `target_resnum` as input
        instead of `target_resname`, if `target_resnum` is given `target_resname`
        is ignored
    target_resname : str, optional
        the residue name of the residue from which you want the nearest
        neighbors, it is suggested to give `target_resnum` as input
        instead of `target_resname`, if `target_resnum` is given `target_resname`
        is ignored
    ignore_resnums : list(int), optional, default=[]
        a list of residue numbers you don't want to include as possible nearest neighbors
    ignore_resnames : list(str), optional,
    default=important_lists.metals + important_lists.trash + ['HOH', 'SOL']
        a list of residue names you don't want to include as possible nearest neighbors,
        NOT case sensitive
    cutoff_angstom : float, optional, default=4.5
        the cutoff to decide if an atom is a neighbor or not, in angstrom
    ignore_hetatms : bool, optional, default=False
        if True will ignore any hetero atom in the PDB (also solvent)

    Returns
    ----------
    NearestNeighborsAtoms
        a class that contains information about the nearest residues' residue names and numbers
        NearestNeighborsAtoms.resnames : list(str)
            the residue names
        NearestNeighborsAtoms.resnumbers : list(int)
            the residue numbers
        the order is consistent

    Raises
    ------------
    ValueError
        if both `target_resname` and `target_resnum` are not given as input
    RuntimeError
        if the target residue is missing

    Notes
    -----------
    If you give `target_resname` and there are more than one residue with this name
    results are unpredictable
    """

    # output class
    class NearestNeighborsAtoms(object):
        """nearest neighbours class
        """
        def __init__(self):

            self.resnames = []

            self.resnumbers = []

        def add(self, resname, resnumber):
            """add a residue
            """

            self.resnames.append(resname)

            self.resnumbers.append(resnumber)

    if ignore_resnums is None:

        ignore_resnums = []

    if ignore_resnames is None:

        ignore_resnames = list(important_lists.metals) +\
             list(important_lists.trash) + ['HOH', 'SOL']

    if target_resnum is not None:

        # to be sure I have a list
        ignore_resnums = list(ignore_resnums)

        ignore_resnums.append(target_resnum)

        check = [target_resnum]

    elif target_resname is not None:

        # in order to deal with possible case issues
        check = [
            target_resname,
            target_resname.upper(),
            target_resname.lower(),
            target_resname.capitalize()
        ]

        # to be sure I have a list
        ignore_resnames = list(ignore_resnames)

        ignore_resnames += check

    else:

        raise ValueError('You must  give a target_resname or a target_resnum')

    #------------------------------------------------------
    #a helper function
    def is_valid_residue(residue):
        """helper function

        needed to check if the given residue shall be checked

        Returns
        ------------
        bool
        """
        valid = True

        #all this things make it false
        if residue.id[1] in ignore_resnums:
            valid = False

        elif residue.resname.strip() in ignore_resnames:
            valid = False

        elif ignore_hetatms:
            if residue.id[0].strip() != '':
                valid = False

        return valid

    #------------------------------------------------------

    # get target residue
    for residue in structure.get_residues():

        if residue.id[1] in check or residue.resname.strip() in check:

            target_residue = residue

            break

    else:

        raise RuntimeError(f'The target residue is missing, {check}')

    nearest_neighbors = NearestNeighborsAtoms()

    #for each residue in the system scan if it has
    #at least 1 atom nearer as the treshold
    #to one atom of the target residue
    for residue in structure.get_residues():

        if is_valid_residue(residue):

            is_neighbor = False

            for atom in residue:

                for target_atom in target_residue:

                    distance = atom.coord - target_atom.coord

                    distance = distance**2

                    distance = np.sum(distance)

                    distance = distance**0.5

                    if distance <= cutoff_angstom:

                        nearest_neighbors.add(resname=residue.resname.strip(),
                                              resnumber=residue.id[1])

                        is_neighbor = True

                        break

                if is_neighbor:

                    break

    return nearest_neighbors


def get_atom_numbers(structure):
    """get a list of atom numbers composing a biopython structure

    it gives you the atom numbers from the atoms composing the given
    biopython structure

    Parameters
    ------------
    structure : bioopython structure or list of atoms

    Returns
    -----------
    list(int)
        the list of all the atom numbers composing the structure

    Notes
    ------------
    useful when definig an area of the structure (like an alchemical region)
    with atom indexes (atom numbers)
    """

    atom_numbers = []

    if hasattr(structure, 'get_atoms'):

        atoms = structure.get_atoms()

    else:  # list of atoms

        atoms = structure

    for atom in atoms:

        atom_numbers.append(atom.get_serial_number())

    return atom_numbers
