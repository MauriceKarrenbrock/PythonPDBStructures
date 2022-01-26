# -*- coding: utf-8 -*-
# pylint: disable=missing-docstring
# pylint: disable=redefined-outer-name
# pylint: disable=no-self-use
# pylint: disable=too-few-public-methods
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

import mdtraj
import numpy as np
from simtk.openmm import unit

import PythonPDBStructures.geometry as geometry
import PythonPDBStructures.pdb.biopython_utils as bio

from . import input_files


class Testget_nearest_neighbors_residues_with_mdtraj():
    def test_works_resnum(self):

        with resources.path(input_files, 'small_dummy.pdb') as pdb:

            with warnings.catch_warnings():

                warnings.simplefilter('ignore')

                traj = mdtraj.load(str(pdb))

        output = geometry.get_nearest_neighbors_residues_with_mdtraj(
            traj, ligand_atoms='resSeq 2', cutoff=0.45 * unit.nanometers)

        assert output == ([0], [1])


class Testget_nearest_neighbors_residues_with_biopython():
    def test_works_resnum(self):

        with resources.path(input_files, 'small_dummy.pdb') as pdb:

            #there are some biopython warnings about the pdb I don't need
            with warnings.catch_warnings():

                warnings.simplefilter('ignore')

                structure = bio.parse_pdb('aaaa', pdb.resolve())

        output = geometry.get_nearest_neighbors_residues_with_biopython(
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


class Testget_inertia_eigenvectors():
    def test_protein_already_in_frame(self):
        with resources.path(input_files,
                            '5rgs_inertia_tensor_frame.pdb') as pdb:

            output_eigenvectors = geometry.get_inertia_eigenvectors(str(pdb))

            # get rotation matrix
            output_eigenvectors = np.linalg.inv(output_eigenvectors)

            # to remove some possible negative values (annoying for testing)
            output_eigenvectors = np.abs(output_eigenvectors)

            expected_eigenvectors = np.array([[1., 0., 0.], [0., 1., 0.],
                                              [0., 0., 1.]])

            assert np.testing.assert_allclose(output_eigenvectors,
                                              expected_eigenvectors,
                                              atol=1.5e-06) is None


class Testget_rotate_coordinates():
    def test_protein_already_in_frame(self):
        with resources.path(input_files,
                            '5rgs_inertia_tensor_frame.pdb') as pdb:

            traj = mdtraj.load(str(pdb))

            identity = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])

            output_xyz = geometry.rotate_coordinates(coordinates=traj.xyz[0],
                                                     rot_matrix=identity,
                                                     check_reflections=False)

            assert np.testing.assert_allclose(output_xyz, traj.xyz[0]) is None
