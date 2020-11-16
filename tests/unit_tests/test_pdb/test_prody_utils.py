# -*- coding: utf-8 -*-
# pylint: disable=missing-docstring
# pylint: disable=redefined-outer-name
# pylint: disable=protected-access
# pylint: disable=no-self-use
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################

from unittest.mock import MagicMock

import PythonPDBStructures.pdb.prody_utils as prody_utils


class Testparse_pdb():
    def test_works(self, mocker):

        m_prody = mocker.patch(
            'PythonPDBStructures.pdb.prody_utils.prody.parsePDB',
            return_value='structure')

        output = prody_utils.parse_pdb('pdb_file.pdb')

        m_prody.assert_called_once_with('pdb_file.pdb')

        assert output == 'structure'


class Testwrite_pdb():
    def test_works(self, mocker):

        m_prody = mocker.patch(
            'PythonPDBStructures.pdb.prody_utils.prody.writePDB')

        prody_utils.write_pdb('TEST_structure', 'pdb_file.pdb')

        m_prody.assert_called_once_with('pdb_file.pdb', 'TEST_structure')


class Testselect():
    def test_works(self):

        structure = MagicMock()

        structure.select.return_value = 'new_structure'

        output = prody_utils.select(structure, 'string')

        structure.select.assert_called_once_with('string')

        assert output == 'new_structure'


class TestProdySelect():
    def test_init(self):

        obj = prody_utils.ProdySelect('TEST_structure')

        assert obj._structure == 'TEST_structure'

    def test_only_protein(self, mocker):

        m_select = mocker.patch('PythonPDBStructures.pdb.prody_utils.select',
                                return_value='NEW_structure')

        obj = prody_utils.ProdySelect('TEST_structure')

        output = obj.only_protein()

        m_select.assert_called_once_with(structure='TEST_structure',
                                         string='protein')

        assert output == 'NEW_structure'

    def test_protein_and_ions_1(self, mocker):

        m_select = mocker.patch('PythonPDBStructures.pdb.prody_utils.select',
                                return_value='NEW_structure')

        obj = prody_utils.ProdySelect('TEST_structure')

        output = obj.protein_and_ions()

        m_select.assert_called_once_with(structure='TEST_structure',
                                         string='protein or ion')

        assert output == 'NEW_structure'

    def test_protein_and_ions_2(self, mocker):

        m_select = mocker.patch('PythonPDBStructures.pdb.prody_utils.select',
                                side_effect=Exception)

        m_only_protein = mocker.patch.object(prody_utils.ProdySelect,
                                             'only_protein',
                                             return_value='NEW_structure')

        obj = prody_utils.ProdySelect('TEST_structure')

        output = obj.protein_and_ions()

        m_select.assert_called_once_with(structure='TEST_structure',
                                         string='protein or ion')

        m_only_protein.assert_called_once()

        assert output == 'NEW_structure'

    def test_resname(self, mocker):

        m_select = mocker.patch('PythonPDBStructures.pdb.prody_utils.select',
                                return_value='NEW_structure')

        obj = prody_utils.ProdySelect('TEST_structure')

        output = obj.resname('test_resname')

        m_select.assert_called_once_with(structure='TEST_structure',
                                         string='resname test_resname')

        assert output == 'NEW_structure'

    def test_resnum(self, mocker):

        m_select = mocker.patch('PythonPDBStructures.pdb.prody_utils.select',
                                return_value='NEW_structure')

        obj = prody_utils.ProdySelect('TEST_structure')

        output = obj.resnum(5)

        m_select.assert_called_once_with(structure='TEST_structure',
                                         string='resnum 5')

        assert output == 'NEW_structure'
