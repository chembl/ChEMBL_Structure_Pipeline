#
#  Copyright (c) 2019 The ChEMBL group
#  All rights reserved.
#
#  This file is part of the ChEMBL_StructurePipeline project.
#  The contents are covered by the terms of the MIT license
#  which is included in the file LICENSE, found at the root
#  of the source tree.
from rdkit import Chem

# Zn not in the list as we have some Zn containing compounds in ChEMBL
# most of them are simple salts
METAL_LIST = [
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Ga", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Cd", "In", "Sn", "La", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Ac",
    "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
    "Yb", "Lu", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
    "Fm", "Md", "No", "Lr", "Ge", "Sb",
]


def exclude_flag(mol, includeRDKitSanitization=True):
    """
    Rules to exclude structures.

    - Metallic or non metallic with more than 7 boron atoms will be excluded
      due to problems when depicting borane compounds.
    """
    rdkit_fails = False
    exclude = False
    metallic = False
    boron_count = 0

    if type(mol) == str:
        mol = Chem.MolFromMolBlock(mol, sanitize=False)
        if includeRDKitSanitization:
            try:
                Chem.SanitizeMol(mol)
            except:
                rdkit_fails = True

    for atom in mol.GetAtoms():
        a_type = atom.GetSymbol()
        if a_type in METAL_LIST:
            metallic = True
        if a_type == "B":
            boron_count += 1

    if metallic or (not metallic and boron_count > 7) or rdkit_fails:
        exclude = True
    return exclude
