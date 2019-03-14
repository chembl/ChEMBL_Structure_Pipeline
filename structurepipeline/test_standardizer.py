#
#  Copyright (c) 2019 Greg Landrum
#  All rights reserved.
#
#  This file is part of the ChEMBL_StructurePipeline project.
#  The contents are covered by the terms of the MIT license
#  which is included in the file LICENSE, found at the root
#  of the source tree.
from . import standardizer
import unittest
from rdkit import Chem


class TestCase(unittest.TestCase):
    def testAtomWiggleBonds(self):
        mb = """
  Mrv1810 03131911492D

  5  4  0  0  0  0            999 V2000
   -5.3125    2.7232    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -4.5980    3.1357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0270    3.1357    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3125    1.8982    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -6.1094    2.5097    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  1  2  1  4  0  0  0
M  END
"""
        outmb = standardizer.standardize_molblock(mb)
        tgtbondblock = """  1  3  1  0
  1  4  1  0
  1  5  1  0
  1  2  1  0"""
        self.assertTrue(outmb.find(tgtbondblock) > 0)

    def testBondWiggleBonds(self):
        mb = """
  Mrv1810 03131911492D          

  4  3  0  0  0  0            999 V2000
   -0.0000    1.8527    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    2.2652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4289    1.8527    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7145    2.2652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  4  1  0  0  0  0
  2  3  1  4  0  0  0
M  END
"""
        outmb = standardizer.standardize_molblock(mb)
        tgtbondblock = """  1  2  2  3
  1  4  1  0
  2  3  1  0
"""
        self.assertTrue(outmb.find(tgtbondblock) > 0)

    def testUncharger1(self):
        mb = """
  Mrv1810 03131914152D          

  4  3  0  0  0  0            999 V2000
   -2.3884   -2.3884    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   -1.6739   -1.9759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9595   -2.3884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2450   -1.9759    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
M  CHG  2   1   1   4  -1
M  END
"""
        m = Chem.MolFromMolBlock(mb, sanitize=False, removeHs=False)
        m.UpdatePropertyCache(False)
        nm = standardizer.uncharge_mol(m)
        print(1, Chem.MolToSmiles(nm))
        for at in nm.GetAtoms():
            self.assertEqual(at.GetFormalCharge(), 0)

    def testUncharger2(self):
        m = Chem.MolFromSmiles("[NH3+]CC[O-]", sanitize=False)
        m.UpdatePropertyCache(False)
        nm = standardizer.uncharge_mol(m)
        print(2, Chem.MolToSmiles(nm))
        for at in nm.GetAtoms():
            self.assertEqual(at.GetFormalCharge(), 0)
