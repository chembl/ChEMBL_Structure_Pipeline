#
#  Copyright (c) 2019 Greg Landrum
#  All rights reserved.
#
#  This file is part of the ChEMBL_StructurePipeline project.
#  The contents are covered by the terms of the MIT license
#  which is included in the file LICENSE, found at the root
#  of the source tree.
import re
from rdkit import Chem
from rdkit import Geometry
from rdkit.Chem import rdinchi
from collections import Counter
from rdkit.Chem.MolStandardize import rdMolStandardize
import rdkit

rdkversion = rdkit.__version__.split(".")[:2]
if rdkversion < ["2019", "03"]:
    raise ValueError("need an RDKit version >= 2019.03.1")


def uncharge_mol(m):
    uncharger = rdMolStandardize.Uncharger()
    return uncharger.uncharge(m)


def standardize_molblock(ctab):
    m = Chem.MolFromMolBlock(ctab, sanitize=False, removeHs=False)

    m = uncharge_mol(m)
    return Chem.MolToMolBlock(m)
