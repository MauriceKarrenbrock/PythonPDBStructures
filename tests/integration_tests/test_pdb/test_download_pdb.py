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

import PythonPDBStructures.pdb.download_pdb as download_pdb


class Testdownload():
    def test_works(self, tmp_path):

        work_dir = tmp_path / 'sub'

        work_dir.mkdir()

        output = download_pdb.download('5aol', 'cif', work_dir)

        assert output == work_dir / '5aol.cif'
