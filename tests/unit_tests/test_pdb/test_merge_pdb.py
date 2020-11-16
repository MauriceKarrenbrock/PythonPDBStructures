# -*- coding: utf-8 -*-
# pylint: disable=missing-docstring
# pylint: disable=redefined-outer-name
# pylint: disable=no-self-use
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################

import PythonPDBStructures.pdb.merge_pdb as merge_pdb


class Testmerge_pdb():
    def test_works(self, mocker):

        pdb_1 = [
            'stuff', 'afafAFGDGSGDSGD',
            'HETATM 3220  C1  UFV A1293      90.951  88.859 -41.631  0.60 14.63           C',
            'ATOM   3221  C1  UFV A1294      90.951  88.859 -41.631  0.60 14.63           C',
            'end stuff'
        ]

        pdb_2 = [
            'more stuff', 'more afafAFGDGSGDSGD',
            'HETATM 3218  NH2 ARG B 290     137.594 109.669 -16.882  1.00 62.21           N',
            'ATOM   3219  NH2 ARG B 291     137.594 109.669 -16.882  1.00 62.21           N',
            'more end stuff'
        ]

        expected = [
            'stuff                                                                         ',
            'afafAFGDGSGDSGD                                                               ',
            'HETATM 3220  C1  UFV A1293      90.951  88.859 -41.631  0.60 14.63           C',
            'ATOM   3221  C1  UFV A1294      90.951  88.859 -41.631  0.60 14.63           C',
            'HETATM 3222  NH2 ARG B1295     137.594 109.669 -16.882  1.00 62.21           N',
            'ATOM   3223  NH2 ARG B1296     137.594 109.669 -16.882  1.00 62.21           N',
            'end stuff                                                                     '
        ]

        m_read = mocker.patch(
            'PythonAuxiliaryFunctions.files_IO.read_file.read_file',
            side_effect=[pdb_1, pdb_2])

        m_write = mocker.patch(
            'PythonAuxiliaryFunctions.files_IO.write_file.write_file')

        merge_pdb.merge_pdb('file_1.pdb', 'file_2.pdb', 'output.pdb')

        m_read.assert_called()

        m_write.assert_called_once_with(lines=expected, file_name='output.pdb')

        ##############
        #PROBLEMA NON AGGIORNO ATOM NUMBER!!!!!!!!
        ##############
