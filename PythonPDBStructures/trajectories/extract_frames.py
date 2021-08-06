# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Functions to extract frames from MD trajectory files
"""

import MDAnalysis as mda


def extract_frames(delta_steps,
                   trajectory,
                   topology,
                   output_name,
                   output_format,
                   starting=0,
                   end=None):
    """Extract frames from a trajectory file at constant `delta_steps`

    Given a MD trajectory file `trajectory` (es .xtc for gromacs) and a
    topology file `topology` (es .tpr for gromacs) it will create a coordinate
    file (pdb, gro, ...) each `delta_steps` MD steps.
    It uses MDAnalysis so check out its documentation for more info
    https://github.com/MDAnalysis/mdanalysis

    Parameters
    ------------
    delta_steps : int
        will print a frame each `delta_steps` MD steps
    trajectory : str or path
        the trajectory file (es xtc for gromacs)
    topology : str or path
        the topology file (es tpr for gromacs)
    output_name : str
        the name of the various coordinates files outputed (pdb, gro,...)
        the complete names will be `output_name`0.`output_format` `output_name`1.`output_format`
        `output_name`2.`output_format` ... `output_name`<n>.`output_format`
    output_format : str
        max 3 letters long it is the file extention of the coordinate file
        you want to get (all the ones that are supported by MDAnalysis are supported)
    starting : int
        the first frame to write, default 0 (the first one)
        useful if you want to jump `starting` initial steps
    end : int
        the last frame to write, default the last one (None)
        useful if you want to jump `end` final steps

    Returns
    ----------
    number_of_created_files : int
        the number of created files, so you know that you will have files with
        names renging from `output_name`0.`output_format` to
        `output_name`(`number_of_created_files` - 1).`output_format`
    """

    #for now MDAnalysis doesn't deal well with pathlib objects
    # issue 2497 n github https://github.com/MDAnalysis/mdanalysis/issues/2497
    topology = str(topology)

    trajectory = str(trajectory)

    universe = mda.Universe(topology, trajectory)

    outputted_file = 0

    for ts in universe.trajectory:

        if end is not None:
            if end < ts.frame:

                break

        if ts.frame < starting:

            continue

        if (ts.frame % delta_steps) == 0:

            universe.atoms.write(
                f'{output_name}{outputted_file}.{output_format}')

            outputted_file += 1

    # I have noticed a memory leack somewhere
    # I haven't been able to localize it yet
    # maybe explicitly cancelling universe may help?
    del universe

    return outputted_file


def extract_all_frames(trajectory, topology, output_name, output_format):
    """Extract ALL frames from a trajectory file

    should be faster than `extract_frames`

    Given a MD trajectory file `trajectory` (es .xtc for gromacs) and a
    topology file `topology` (es .tpr for gromacs) it will create a coordinate
    file (pdb, gro, ...) for each single frame.
    It uses MDAnalysis so check out its documentation for more info
    https://github.com/MDAnalysis/mdanalysis

    Parameters
    ------------
    trajectory : str or path
        the trajectory file (es xtc for gromacs)
    topology : str or path
        the topology file (es tpr for gromacs)
    output_name : str
        the name of the various coordinates files outputed (pdb, gro,...)
        the complete names will be `output_name`0.`output_format` `output_name`1.`output_format`
        `output_name`2.`output_format` ... `output_name`<n>.`output_format`
    output_format : str
        max 3 letters long it is the file extention of the coordinate file
        you want to get (all the ones that are supported by MDAnalysis are supported)

    Returns
    ----------
    number_of_created_files : int
        the number of created files, so you know that you will have files with
        names renging from `output_name`0.`output_format` to
        `output_name`(`number_of_created_files` - 1).`output_format`
    """

    #for now MDAnalysis doesn't deal well with pathlib objects
    # issue 2497 n github https://github.com/MDAnalysis/mdanalysis/issues/2497
    topology = str(topology)

    trajectory = str(trajectory)

    universe = mda.Universe(topology, trajectory)

    for ts in universe.trajectory:

        universe.atoms.write(f'{output_name}{ts.frame}.{output_format}')

    n_frames = universe.trajectory.n_frames

    # I have noticed a memory leack somewhere
    # I haven't been able to localize it yet
    # maybe explicitly cancelling universe may help?
    del universe

    return n_frames
