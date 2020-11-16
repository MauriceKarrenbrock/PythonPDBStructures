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

from pathlib import Path

import Bio.PDB.PDBList
import pytest

import PythonPDBStructures.pdb.download_pdb as download_pdb


class Testdownload():
    def test_raises_ValueError(self):

        with pytest.raises(ValueError):

            download_pdb.download('5aol', 'WRONG FILE TYPE')

    def test_works_cif(self, mocker):

        m_retrieve = mocker.patch.object(Bio.PDB.PDBList,
                                         'retrieve_pdb_file',
                                         return_value='TEST_PATH')

        m_exists = mocker.patch.object(Path, 'exists', return_value=True)

        m_cwd = mocker.patch.object(Path, 'cwd', return_value='WORK_PATH')

        output = download_pdb.download('5aol', 'cif')

        m_retrieve.assert_called_once_with('5aol',
                                           False,
                                           'WORK_PATH',
                                           file_format='mmCif',
                                           overwrite=True)

        m_exists.assert_called_once()

        m_cwd.assert_called_once()

        assert isinstance(output, Path)

    def test_works_pdb(self, mocker):

        m_retrieve = mocker.patch.object(Bio.PDB.PDBList,
                                         'retrieve_pdb_file',
                                         return_value='TEST_PATH')

        m_exists = mocker.patch.object(Path, 'exists', return_value=True)

        m_cwd = mocker.patch.object(Path, 'cwd')

        output = download_pdb.download('5aol', 'pdb', 'SOME_WORK_DIR')

        m_retrieve.assert_called_once_with('5aol',
                                           False,
                                           'SOME_WORK_DIR',
                                           file_format='pdb',
                                           overwrite=True)

        m_exists.assert_called_once()

        m_cwd.assert_not_called()

        assert isinstance(output, Path)

    def test_raises_FileNotFoundError(self, mocker):

        m_retrieve = mocker.patch.object(Bio.PDB.PDBList,
                                         'retrieve_pdb_file',
                                         return_value='TEST_PATH')

        m_exists = mocker.patch.object(Path, 'exists', return_value=False)

        m_cwd = mocker.patch.object(Path, 'cwd')

        with pytest.raises(FileNotFoundError):

            download_pdb.download('5aol', 'pdb', 'SOME_WORK_DIR')

        m_retrieve.assert_called_once_with('5aol',
                                           False,
                                           'SOME_WORK_DIR',
                                           file_format='pdb',
                                           overwrite=True)

        m_exists.assert_called_once()

        m_cwd.assert_not_called()
