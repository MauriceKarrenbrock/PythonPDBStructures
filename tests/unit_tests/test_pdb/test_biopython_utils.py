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

import unittest.mock as mock

import Bio.PDB
import pytest

import PythonPDBStructures.pdb.biopython_utils as bio


class Testparse_pdb():
    def test_works(self, mocker):

        prot_id = '5aol'

        file_name = 'name.pdb'

        m_parser = mocker.patch.object(Bio.PDB.PDBParser, 'get_structure')

        output = bio.parse_pdb(prot_id, file_name)

        m_parser.assert_called_once_with(prot_id, file_name)

        assert isinstance(output, mock.MagicMock)


class Testparse_mmcif():
    def test_works(self, mocker):

        prot_id = '5aol'

        file_name = 'name.cif'

        m_parser = mocker.patch.object(Bio.PDB.MMCIFParser, 'get_structure')

        output = bio.parse_mmcif(prot_id, file_name)

        m_parser.assert_called_once_with(prot_id, file_name)

        assert isinstance(output, mock.MagicMock)


class Testmmcif2dict():
    def test_works(self, mocker):

        file_name = 'name.cif'

        m_parser = \
        mocker.patch('PythonPDBStructures.pdb.biopython_utils.Bio.PDB.MMCIF2Dict.MMCIF2Dict')

        output = bio.mmcif2dict(file_name)

        m_parser.assert_called_once_with(file_name)

        assert isinstance(output, mock.MagicMock)


class Testwrite_pdb():
    def test_works_is_instance(self, mocker):

        m_isinstance = \
        mocker.patch('PythonPDBStructures.pdb.biopython_utils.isinstance', return_value=True)

        m_hasattr = \
        mocker.patch('PythonPDBStructures.pdb.biopython_utils.hasattr')

        m_struct = mocker.patch.object(Bio.PDB.PDBIO, 'set_structure')

        m_save = mocker.patch.object(Bio.PDB.PDBIO, 'save')

        bio.write_pdb('fake_structure_object', 'file_name.pdb')

        m_isinstance.assert_called_once()

        m_hasattr.assert_not_called()

        m_struct.assert_called_once_with('fake_structure_object')

        m_save.assert_called_once_with('file_name.pdb')

    def test_works_hasattr(self, mocker):
        class FakeAtom(object):
            def __init__(self):

                self.level = 'A'

        fake_atom_1 = FakeAtom()
        fake_atom_2 = FakeAtom()

        fake_structure = (fake_atom_1, fake_atom_2)

        m_isinstance = \
        mocker.patch('PythonPDBStructures.pdb.biopython_utils.isinstance', return_value=False)

        m_hasattr = \
        mocker.patch('PythonPDBStructures.pdb.biopython_utils.hasattr', return_value=True)

        m_struct = mocker.patch.object(Bio.PDB.PDBIO, 'set_structure')

        m_save = mocker.patch.object(Bio.PDB.PDBIO, 'save')

        bio.write_pdb(fake_structure, 'file_name.pdb')

        m_isinstance.assert_called_once()

        m_hasattr.assert_called_once()

        m_struct.assert_called_once_with(fake_structure)

        m_save.assert_called_once_with('file_name.pdb')

    def test_raises(self, mocker):
        class FakeAtom(object):
            def __init__(self):

                self.level = 'A'

        fake_atom_1 = FakeAtom()
        fake_atom_2 = FakeAtom()

        fake_structure = (fake_atom_1, fake_atom_2)

        m_isinstance = \
        mocker.patch('PythonPDBStructures.pdb.biopython_utils.isinstance', return_value=False)

        m_hasattr = \
        mocker.patch('PythonPDBStructures.pdb.biopython_utils.hasattr', return_value=False)

        with pytest.raises(TypeError):

            bio.write_pdb(fake_structure, 'file_name.pdb')

        m_isinstance.assert_called_once()

        m_hasattr.assert_called_once()


class Testwrite_mmcif():
    def test_works_is_instance(self, mocker):

        m_isinstance = \
        mocker.patch('PythonPDBStructures.pdb.biopython_utils.isinstance', return_value=True)

        m_hasattr = \
        mocker.patch('PythonPDBStructures.pdb.biopython_utils.hasattr')

        m_struct = mocker.patch.object(Bio.PDB.MMCIFIO, 'set_structure')

        m_save = mocker.patch.object(Bio.PDB.MMCIFIO, 'save')

        bio.write_mmcif('fake_structure_object', 'file_name.cif')

        m_isinstance.assert_called_once()

        m_hasattr.assert_not_called()

        m_struct.assert_called_once_with('fake_structure_object')

        m_save.assert_called_once_with('file_name.cif')

    def test_works_hasattr(self, mocker):
        class FakeAtom(object):
            def __init__(self):

                self.level = 'A'

        fake_atom_1 = FakeAtom()
        fake_atom_2 = FakeAtom()

        fake_structure = (fake_atom_1, fake_atom_2)

        m_isinstance = \
        mocker.patch('PythonPDBStructures.pdb.biopython_utils.isinstance', return_value=False)

        m_hasattr = \
        mocker.patch('PythonPDBStructures.pdb.biopython_utils.hasattr', return_value=True)

        m_struct = mocker.patch.object(Bio.PDB.MMCIFIO, 'set_structure')

        m_save = mocker.patch.object(Bio.PDB.MMCIFIO, 'save')

        bio.write_mmcif(fake_structure, 'file_name.cif')

        m_isinstance.assert_called_once()

        m_hasattr.assert_called_once()

        m_struct.assert_called_once_with(fake_structure)

        m_save.assert_called_once_with('file_name.cif')

    def test_raises(self, mocker):
        class FakeAtom(object):
            def __init__(self):

                self.level = 'A'

        fake_atom_1 = FakeAtom()
        fake_atom_2 = FakeAtom()

        fake_structure = (fake_atom_1, fake_atom_2)

        m_isinstance = \
        mocker.patch('PythonPDBStructures.pdb.biopython_utils.isinstance', return_value=False)

        m_hasattr = \
        mocker.patch('PythonPDBStructures.pdb.biopython_utils.hasattr', return_value=False)

        with pytest.raises(TypeError):

            bio.write_mmcif(fake_structure, 'file_name.cif')

        m_isinstance.assert_called_once()

        m_hasattr.assert_called_once()


class Testwrite_dict2mmcif():
    def test_works(self, mocker):

        fake_dict = {}

        m_dict = mocker.patch.object(Bio.PDB.MMCIFIO, 'set_dict')

        m_save = mocker.patch.object(Bio.PDB.MMCIFIO, 'save')

        bio.write_dict2mmcif(fake_dict, 'file_name.cif')

        m_dict.assert_called_once_with(fake_dict)

        m_save.assert_called_once_with('file_name.cif')
