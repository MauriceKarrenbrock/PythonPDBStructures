# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Contains the functions to parse and write (PDB MMCIF) files with biopython
"""

import Bio.PDB
import Bio.PDB.Entity
import Bio.PDB.MMCIF2Dict

######################################
# Parsers
#####################################


def parse_pdb(protein_id, file_name):
    """Parse a PDB file with Biopython

    This function uses Biopython Bio.PDB.PDBParser to
    return a Biopython structure from a PDB file.

    Parameters
    -------------
    protein_id : str
        the PDB databank id
    file_name : str
        file name (or path) of the PDB file

    Returns
    --------------
    structure
        Biopython structure
    """

    parser = Bio.PDB.PDBParser()

    structure = parser.get_structure(protein_id, file_name)

    return structure


def parse_mmcif(protein_id, file_name):
    """Parse a PDB file with Biopython

    This function uses Biopython Bio.PDB.MMCIFParser to
    return a Biopython structure from a PDB file.

    Parameters
    -------------
    protein_id : str
        the PDB databank id
    file_name : str
        file name (or path) of the MMCif file

    Returns
    --------------
    structure
        Biopython structure
    """

    parser = Bio.PDB.MMCIFParser()

    structure = parser.get_structure(protein_id, file_name)

    return structure


def mmcif2dict(file_name):
    """Uses Bio.PDB.MMCIF2DICT to return a dictionary

    Uses Bio.PDB.MMCIF2DICT to return a dictionary of the mmcif file

    Parameters
    -------------
    file_name : str
        file name (or path) of the MMCif file

    Returns
    -------------------
    dict
    """

    return Bio.PDB.MMCIF2Dict.MMCIF2Dict(file_name)


##################################################
# Writers
##################################################


def write_pdb(structure, file_name='file.pdb'):
    """writes a PDB file when given a Biopython structure

    Parameters
    ------------
    structure : Bio.PDB.Entity.Entity
    file_name :: str
        default file.pdb
    """

    # Structure, Model, Chain, Residue
    if isinstance(structure, Bio.PDB.Entity.Entity):
        pass
    # List of Atoms
    elif hasattr(structure,
                 '__iter__') and [x for x in structure if x.level == 'A']:
        pass
    else:  # Some other weirdo object
        raise TypeError('Need a Bio.PDB.Entity instance like:\n'
                        'Structure, Model, Chain, Residue, list of Atoms.')

    writer = Bio.PDB.PDBIO()

    writer.set_structure(structure)

    writer.save(file_name)


def write_mmcif(structure, file_name='file.cif'):
    """writes a MMCIF file when given a Biopython structure

    Parameters
    ------------
    structure : Bio.PDB.Entity.Entity
    file_name :: str
        default file.pdb
    """
    # Structure, Model, Chain, Residue
    if isinstance(structure, Bio.PDB.Entity.Entity):
        pass
    # List of Atoms
    elif hasattr(structure,
                 '__iter__') and [x for x in structure if x.level == 'A']:
        pass
    else:  # Some other weirdo object
        raise TypeError('Need a Bio.PDB.Entity instance like:\n'
                        'Structure, Model, Chain, Residue, list of Atoms.')

    writer = Bio.PDB.MMCIFIO()

    writer.set_structure(structure)

    writer.save(file_name)


def write_dict2mmcif(dictionary, file_name='file.cif'):
    """Writes a mmcif file starting from a dictionary

    Writes a mmcif file starting from a dictionary
    obtained from mmcif2dict (that uses Bio.PDB.MMCIF2DICT )

    Parameters
    ---------------
    dictionary : dict
        a dictionary containing all the mmcif infos, obtained with mmcif2dict
    file_name : str
    """

    p = Bio.PDB.MMCIFIO()
    p.set_dict(dictionary)
    p.save(file_name)
