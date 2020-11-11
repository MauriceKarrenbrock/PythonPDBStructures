# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Add chain ID
"""

import PythonAuxiliaryFunctions.files_IO.read_file as read_file
import PythonAuxiliaryFunctions.files_IO.write_file as write_file


def add_chain_id_pdb(pdb_file, chain='A'):
    """Adds a chain ID to a PDB file

    This is a patch because some MD programs remove the chain id from
    pdb files and this confuses some pdb parsers (works on PDB files only)

    Parameters
    ---------------
    pdb_file : str
        the pdb file to edit (name or path)
    chain : str
        default 'A', the chain id to add to the pdb_file
    """

    chain = chain.upper().strip()

    lines = read_file.read_file(file_name=pdb_file)

    for i in range(len(lines)):

        if lines[i][0:4] == 'ATOM' or lines[i][0:6] == 'HETATM' or lines[i][
                0:3] == 'TER':

            lines[i] = lines[i][:20] + '{0:>2}'.format(chain) + lines[i][22:]

            lines[i] = lines[i].strip('\n') + '\n'

    write_file.write_file(lines=lines, file_name=pdb_file)
