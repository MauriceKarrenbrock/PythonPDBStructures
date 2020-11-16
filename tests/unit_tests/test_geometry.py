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

import PythonPDBStructures.geometry as geometry
import PythonPDBStructures.important_lists


class Testget_center_of_mass():
    @pytest.mark.parametrize('test_type, geometric, expected_COM',
                             [('geometric_false', False, (-1. / 3., 0., 0.)),
                              ('geometric_true', True, (0., 0., 0.))])
    def test_atom_list(self, mocker, test_type, geometric, expected_COM):

        print('Logging test type for visibility: ' + test_type)

        test_weights = {'C': 2. / 3., 'He': 1. / 3.}

        class Atom(object):
            def __init__(self, element, coord):

                self.element = element
                self.coord = coord

        atom_1 = Atom('c', (-1., 0., 0.))

        atom_2 = Atom('he', (1., 0., 0.))

        molecule = (atom_1, atom_2)

        mocker.patch.dict(PythonPDBStructures.important_lists.atom_weights,
                          test_weights,
                          clear=True)

        output_COM = geometry.get_center_of_mass(molecule, geometric)

        assert output_COM == expected_COM

    def test_raises(self, mocker):

        test_weights = {}

        class Atom(object):
            def __init__(self, element, coord):

                self.element = element
                self.coord = coord

        atom_1 = Atom('c', (-1., 0., 0.))

        atom_2 = Atom('he', (1., 0., 0.))

        molecule = (atom_1, atom_2)

        mocker.patch.dict(PythonPDBStructures.important_lists.atom_weights,
                          test_weights,
                          clear=True)

        with pytest.raises(ValueError):
            geometry.get_center_of_mass(molecule, False)

    def test_biopython_entity(self, mocker):

        test_weights = {'C': 2. / 3., 'He': 1. / 3.}

        class Atom(object):
            def __init__(self, element, coord):

                self.element = element
                self.coord = coord

        atom_1 = Atom('c', (-1., 0., 0.))

        atom_2 = Atom('he', (1., 0., 0.))

        class Molecule(object):
            def __init__(self, atoms):

                self.atoms = atoms

            def get_atoms(self):

                return self.atoms

        molecule = Molecule((atom_1, atom_2))

        mocker.patch('PythonPDBStructures.geometry.isinstance',
                     return_value=True)

        mocker.patch.dict(PythonPDBStructures.important_lists.atom_weights,
                          test_weights,
                          clear=True)

        output_COM = geometry.get_center_of_mass(molecule, False)

        assert output_COM == (-1. / 3., 0., 0.)