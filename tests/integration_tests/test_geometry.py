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

try:
    import importlib.resources as resources  # python > 3.6
except ImportError:
    import importlib_resources as resources

import warnings

import PythonPDBStructures.geometry as geometry
import PythonPDBStructures.pdb.biopython_utils as bio

from . import input_files


class Testget_nearest_neighbors_residues():
    def test_works_resnum(self):

        with resources.path(input_files, 'small_dummy.pdb') as pdb:

            #there are some biopython warnings about the pdb I don't need
            with warnings.catch_warnings():

                warnings.simplefilter('ignore')

                structure = bio.parse_pdb('aaaa', pdb.resolve())

        output = geometry.get_nearest_neighbors_residues(
            structure,
            target_resnum=1,
            target_resname=None,
            ignore_resnums=None,
            ignore_resnames=['HOH'],
            cutoff_angstom=4.5,
            ignore_hetatms=False)

        assert output.resnames == ['B2']

        assert output.resnumbers == [2]


class Testget_atom_numbers():
    def test_works(self):

        with resources.path(input_files, 'small_dummy.pdb') as pdb:

            #there are some biopython warnings about the pdb I don't need
            with warnings.catch_warnings():

                warnings.simplefilter('ignore')

                structure = bio.parse_pdb('aaaa', pdb.resolve())

        expected_output = list(range(1, 4492))

        output = geometry.get_atom_numbers(structure)

        assert output == expected_output
