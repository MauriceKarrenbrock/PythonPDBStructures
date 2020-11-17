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

from pathlib import Path

import PythonAuxiliaryFunctions.files_IO.read_file as read_file

import PythonPDBStructures.trajectories.extract_frames as extract_frames


class Testextract_frames():
    def test_works(self, tmp_path):

        input_dir = Path('tests/integration_tests/input_files')

        trajectory = input_dir / 'ch2cl2.trr'

        topology = input_dir / 'ch2cl2.tpr'

        output_path = tmp_path / 'ch2cl2_extract_frames'

        output_path.mkdir()

        output_file = output_path / 'output'

        output = extract_frames.extract_frames(10, trajectory, topology,
                                               output_file, 'pdb', 3, 89)

        assert output == 8

        assert set(output_path.iterdir()) == set([
            output_path / 'output0.pdb', output_path / 'output1.pdb',
            output_path / 'output2.pdb', output_path / 'output3.pdb',
            output_path / 'output4.pdb', output_path / 'output5.pdb',
            output_path / 'output6.pdb', output_path / 'output7.pdb'
        ])

        assert read_file.read_file(output_path / 'output0.pdb')[1:] != \
            read_file.read_file(output_path / 'output1.pdb')[1:]
