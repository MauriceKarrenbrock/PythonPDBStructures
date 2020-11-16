# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Functions to download the pdb and mmcif files from the wwPDB
"""

from pathlib import Path

import Bio.PDB


def download(protein_id, file_type='cif', pdir=None):
    """The function downloads a PDB or a mmCIF from wwwPDB in a selected directory

    the default directory is the working directory
    it returns the filename (str)

    Parameters
    --------------
    protein_id : str
        it is the protein to download
    file_type : str
        it can be pdb or cif depending on the format required, default cif
    pdir : str
        default is working directory
        it is the directory where the file is saved

    Returns
    ----------
    file_path : pathlib.Path
        check pathlib documentation for more info on this kind of
        path object

    Raises
    ------------
    FileNotFoundError
        if the file is not downloaded correctly
    """

    if pdir is None:

        pdir = Path.cwd()

    if file_type == 'cif':
        file_type = 'mmCif'

    elif file_type not in ('pdb', 'cif', 'mmCif'):
        raise ValueError(f"Must be 'pdb' or 'cif' not {file_type}")

    pdbl = Bio.PDB.PDBList()
    file_path = pdbl.retrieve_pdb_file(protein_id,
                                       False,
                                       pdir,
                                       file_format=file_type,
                                       overwrite=True)

    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(
            f'Was not able to download the protein or to find {file_path}')

    return file_path
