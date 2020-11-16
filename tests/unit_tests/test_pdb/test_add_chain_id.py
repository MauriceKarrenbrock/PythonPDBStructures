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

import pytest

import PythonPDBStructures.pdb.add_chain_id as add_chain_id


class Testadd_chain_id():
    @pytest.mark.parametrize('test_type, chain', [('1_letter_chain', 'b'),
                                                  ('2_letter_chain', 'bb')])
    def test_works(self, mocker, test_type, chain):

        print('Logging test type for visibility: ' + test_type)

        pdb_file = [
            'AAAAAAAAAA\n',
            'ATOM    698  C   VAL    91      15.434   8.943  21.950  1.00 16.60           C\n',
            'HETATM  698  C   VAL    91      15.434   8.943  21.950  1.00 16.60           C\n',
            'TER     698  C   VAL    91      15.434   8.943  21.950  1.00 16.60           C\n',
            'GSGFGSGSGASG'
        ]

        expected = [
            'AAAAAAAAAA\n',
            'ATOM    698  C   VAL{:>2}  91      15.434   8.943  21.950  1.00 16.60           C\n'
            .format(chain.upper()),  # pylint: disable=line-too-long
            'HETATM  698  C   VAL{:>2}  91      15.434   8.943  21.950  1.00 16.60           C\n'
            .format(chain.upper()),  # pylint: disable=line-too-long
            'TER     698  C   VAL{:>2}  91      15.434   8.943  21.950  1.00 16.60           C\n'
            .format(chain.upper()),  # pylint: disable=line-too-long
            'GSGFGSGSGASG'
        ]

        m_read = mocker.patch(
            'PythonAuxiliaryFunctions.files_IO.read_file.read_file',
            return_value=pdb_file)

        m_write = mocker.patch(
            'PythonAuxiliaryFunctions.files_IO.write_file.write_file')

        add_chain_id.add_chain_id_pdb('file.pdb', chain)

        m_read.assert_called_once()

        m_write.assert_called_once_with(lines=expected, file_name='file.pdb')

    def test_raises(self):

        with pytest.raises(ValueError):

            add_chain_id.add_chain_id_pdb('file.pdb', 'AAA')
