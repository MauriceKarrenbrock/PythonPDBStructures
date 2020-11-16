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

import unittest.mock as mock
from unittest.mock import MagicMock

import PythonPDBStructures.trajectories.extract_frames as extract_frames


class Testextract_frames():
    def test_works(self, mocker):
        class Frame():
            def __init__(self, frame):

                self.frame = frame

        #create test trajectory
        test_trajectory = []
        for i in range(100):

            test_trajectory.append(Frame(i))

        test_trajectory = tuple(test_trajectory)

        universe = mocker.patch('MDAnalysis.Universe')

        universe.return_value = universe

        atoms = MagicMock()

        universe.select_atoms.return_value = atoms

        universe.trajectory = test_trajectory

        output = extract_frames.extract_frames(10, 'trajectory', 'topology',
                                               'output', 'pdb', 3, 89)

        assert output == 8

        assert atoms.write.call_args_list == [
            mock.call('output0.pdb'),
            mock.call('output1.pdb'),
            mock.call('output2.pdb'),
            mock.call('output3.pdb'),
            mock.call('output4.pdb'),
            mock.call('output5.pdb'),
            mock.call('output6.pdb'),
            mock.call('output7.pdb')
        ]

        universe.assert_called()
