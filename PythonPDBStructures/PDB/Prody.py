# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Contains the functions to write and parse pdb files with prody
"""

import prody


def parse_pdb(file_name):
    """Parses a PDB file with ProDy and returns a ProDy structure (prody.AtomGroup)

    Parameters
    -----------
    file_name : str

    Returns
    ----------
    structure : prody.AtomGroup
    """

    structure = prody.parsePDB(file_name)

    return structure


def write_pdb(structure, file_name='file.pdb'):
    """Takes a Prody structure prody.AtomGroup and writes it on a pdb file

    Parameters
    -------------
    strucure : prody.AtomGroup
    file_name : str
        default "file.pdb"
    """

    prody.writePDB(file_name, structure)


def select(structure, string):
    """Uses the Prody select function with string as command

    Makes you easily extract some parts of a structure

    Parameters
    --------------
    structure : prody.AtomGroup
    string : str
        this is the command that will be passed to prody select

    Returns
    ---------
    new_structure : prody.AtomGroup
    """

    new_structure = structure.select(string)

    return new_structure


class ProdySelect(object):
    """This cass is a smart facade for the select function

    Parameters
    -----------
    structure : prody.AtomGroup
    """
    def __init__(self, structure):

        self._structure = structure

    def only_protein(self):
        """Returns a prody structure containing only the protein (no ions)

        Returns
        -----------
        prody.AtomGroup
        """

        return select(structure=self._structure, string='protein')

    def protein_and_ions(self):
        """Returns a prody structure containing only the protein and the inorganic ions (if any)

        Returns
        -----------
        prody.AtomGroup
        """

        try:

            return select(structure=self._structure, string='protein or ion')

        # pylint: disable=bare-except
        except:

            return self.only_protein()

    def resname(self, resname):
        """returns a Prody structure only containing any residue with that residue name

        Parameters
        ------------
        resname : str
            residue name

        Returns
        -----------
        prody.AtomGroup
        """

        resname = resname.strip().upper()

        return select(structure=self._structure, string='resname ' + resname)

    def resnum(self, resnum):
        """returns a Prody structure only containing the residue with that residue number

        Parameters
        -----------
        resnum : int
            the residue number

        Returns
        -----------
        prody.AtomGroup
        """

        return select(structure=self._structure, string=f'resnum {resnum}')
