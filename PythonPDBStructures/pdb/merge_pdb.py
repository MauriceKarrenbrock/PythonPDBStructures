# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""
This file contains the functions to merge PDB or mmCIF
They are useful wen you need to merge one or more organic ligands with a protein
"""

from PythonAuxiliaryFunctions.files_IO import read_file, write_file


def merge_pdb(input_pdb_1, input_pdb_2, output_pdb):
    """Merge 2 PDB files

    this function is brutal and memory consuming I should do it better in the future
    it adds file 2 to file 1 (the output file will be in order 1-2)

    The header of file 2 will be omitted and that of file 1 untouched

    Parameters
    ------------
    input_pdb_1 : str
    input_pdb_2 : str
    output_pdb : str
        can also be one of the input ones
    """

    lines_1 = read_file.read_file(input_pdb_1)

    #get the index of the line with the last ATOM HETATM or TER line
    #and get the resnum of this last residue
    for i in range(len(lines_1) - 1, -1, -1):

        # pylint: disable=no-else-break
        if lines_1[i][0:4] == 'ATOM' or lines_1[i][0:6] == 'HETATM':

            residue_number = int(lines_1[i][22:26].strip())

            index_protein_file = i + 1

            break

        elif lines_1[i][0:3] == 'TER':

            #some TER lines are non standard and don't contain the residue number
            residue_number = int(lines_1[i - 1][22:26].strip())

            index_protein_file = i + 1

            break

    lines_2 = read_file.read_file(input_pdb_2)
    #find first coordinate line in input_pdb_2
    for i in range(len(lines_2)):

        if lines_2[i][0:4] == 'ATOM' or lines_2[i][0:6] == 'HETATM':

            beginnig_pdb_2 = i

            break

    for i in range(len(lines_2), -1, -1):

        if lines_2[i][0:4] == 'ATOM' or lines_2[i][0:6] == 'HETATM' or lines_2[
                i][0:3] == 'TER':

            end_pdb_2 = i

            break

    lines_2 = lines_2[beginnig_pdb_2:end_pdb_2 + 1]

    residue_number += 1

    for i in range(len(lines_2)):

        lines_2[i] = lines_2[i][:22] + '{0:>4}'.format(
            residue_number) + lines_2[i][26:]

    #insert the ligands in the right place of the protein_file list
    lines_1[index_protein_file:index_protein_file] = lines_2

    write_file.write_file(lines=lines_1, file_name=output_pdb)