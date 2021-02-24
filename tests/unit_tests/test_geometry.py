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

import numpy as np
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
            def __init__(self, element, coord, name='DUM'):

                self.element = element
                self.coord = coord
                self.name = name

        atom_1 = Atom('c', (-1., 0., 0.))

        atom_2 = Atom('he', (1., 0., 0.))

        molecule = (atom_1, atom_2)

        mocker.patch.dict(PythonPDBStructures.important_lists.atom_weights,
                          test_weights,
                          clear=True)

        output_COM = geometry.get_center_of_mass(molecule, geometric)

        assert np.testing.assert_allclose(output_COM,
                                          np.array(expected_COM)) is None

    def test_raises(self, mocker):

        test_weights = {}

        class Atom(object):
            def __init__(self, element, coord, name='DUM'):

                self.element = element
                self.coord = coord
                self.name = name

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
            def __init__(self, element, coord, name='DUM'):

                self.element = element
                self.coord = coord
                self.name = name

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

        assert np.testing.assert_allclose(output_COM,
                                          np.array((-1. / 3., 0., 0.))) is None

    def test_backup_name(self):
        class Atom(object):
            def __init__(self, element, coord, name):

                self.element = element
                self.coord = coord
                self.name = name

        atom_1 = Atom('Does_not_Exist', (-1., 0., 0.),
                      '13132132beeeeeeerfferferf3235525')  #Be

        atom_2 = Atom('Does_not_Exist', (1., 0., 0.),
                      '13132132beeeeeeerfferferf3235525')  #Be

        molecule = (atom_1, atom_2)

        output_COM = geometry.get_center_of_mass(molecule, False)

        assert np.testing.assert_allclose(output_COM, np.array(
            (0., 0., 0.))) is None


class Testget_nearest_neighbors_residues():
    def test_works_resnum(self):
        class Atom(object):
            def __init__(self, element, coord, name='atom'):

                self.element = element
                self.coord = coord
                self.name = name

        class Residue(object):
            def __init__(self, atom_list, resnum, resname='RES'):

                self.atom_list = atom_list

                self.id = ['', resnum]

                self.resname = resname

            def get_atoms(self):

                return self.atom_list

            def __iter__(self):

                return self.atom_list.__iter__()

        class Structure(object):
            def __init__(self, residue_list):

                self.residue_list = residue_list

            def get_residues(self):

                return self.residue_list

        atom_1 = Atom('c', np.array([-1., 0., 0.]))

        atom_2 = Atom('he', np.array([1., 0., 0.]))

        atom_3 = Atom('he', np.array([10., 10., 10.]))

        residue_near = Residue([atom_1, atom_3], 1, '1')

        residue_far = Residue([atom_3, atom_3], 2, '2')

        residue_ignore = Residue([atom_1, atom_3], 3, 'NO')

        residue_target = Residue([atom_2, atom_2], 4, '4')

        structure = Structure(
            [residue_near, residue_far, residue_ignore, residue_target])

        output = geometry.get_nearest_neighbors_residues(
            structure,
            target_resnum=4,
            target_resname=None,
            ignore_resnums=None,
            ignore_resnames=['NO'],
            cutoff_angstom=4.5,
            ignore_hetatms=False)

        assert output.resnames == ['1']

        assert output.resnumbers == [1]


class Testget_atom_numbers():
    def test_works(self):
        class Atom(object):
            def __init__(self, element, coord, name='atom', seria_number=1):

                self.element = element
                self.coord = coord
                self.name = name

                self.seria_number = seria_number

            def get_serial_number(self):

                return self.seria_number

        class Residue(object):
            def __init__(self, atom_list, resnum, resname='RES'):

                self.atom_list = atom_list

                self.id = ['', resnum]

                self.resname = resname

            def get_atoms(self):

                return self.atom_list

            def __iter__(self):

                return self.atom_list.__iter__()

        class Structure(object):
            def __init__(self, residue_list):

                self.residue_list = residue_list

            def get_residues(self):

                return self.residue_list

            def get_atoms(self):

                atoms = []

                for i in self.residue_list:

                    atoms += i.get_atoms()

                return atoms

        atom_1 = Atom('c', np.array([-1., 0., 0.]), seria_number=1)

        atom_2 = Atom('he', np.array([1., 0., 0.]), seria_number=2)

        atom_3 = Atom('he', np.array([10., 10., 10.]), seria_number=3)

        residue_1 = Residue([atom_1, atom_3], 1, '1')

        residue_2 = Residue([atom_3, atom_3], 2, '2')

        residue_3 = Residue([atom_1, atom_3], 3, 'NO')

        residue_4 = Residue([atom_2, atom_2], 4, '4')

        structure = Structure([residue_1, residue_2, residue_3, residue_4])

        expected_output = [1, 3, 3, 3, 1, 3, 2, 2]

        output = geometry.get_atom_numbers(structure)

        assert output == expected_output
