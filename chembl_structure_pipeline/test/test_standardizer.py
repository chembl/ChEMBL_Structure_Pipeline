#
#  Copyright (c) 2019 Greg Landrum
#  All rights reserved.
#
#  This file is part of the ChEMBL_StructurePipeline project.
#  The contents are covered by the terms of the MIT license
#  which is included in the file LICENSE, found at the root
#  of the source tree.
from .. import standardizer
import unittest
import math
import os
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdqueries
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.info")
_test_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")

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
        for at in nm.GetAtoms():
            self.assertEqual(at.GetFormalCharge(), 0)

        m = Chem.MolFromSmiles("[NH3+]CC[O-]", sanitize=False)
        m.UpdatePropertyCache(False)
        nm = standardizer.uncharge_mol(m)
        for at in nm.GetAtoms():
            self.assertEqual(at.GetFormalCharge(), 0)

    def testUncharger2(self):
        data = (
            ('[NH3+]CC[O-]', 'NCCO'),
            ('[NH3+]CCC[O-].[Na+]', 'NCCC[O-].[Na+]'),
            ('[Na+].[NH3+]CCC[O-]', 'NCCC[O-].[Na+]'),
            ('[NH3+]CCO', 'NCCO'),
            ('NCC[O-]', 'NCCO'),
            ('[Cl-].[NH3+]CCC[O-]', 'Cl.NCCCO'),
            ('[N+](C)(C)(C)CCC[O-]', '[N+](C)(C)(C)CCC[O-]'),
            ('[NH3+]CC([O-])C[O-]', 'NCC(O)CO'),
            # ('[NH3+]CC([O-])C[O-].[Na+]','NCC(O)C[O-].[Na+]'),
            ('[NH3+]CCC[O-].[NH+](C)(C)C', 'CN(C)C.NCCCO'))
        for ismi, esmi in data:
            esmi = Chem.CanonSmiles(esmi)
            m = Chem.MolFromSmiles(ismi, sanitize=False)
            m.UpdatePropertyCache(False)
            nm = standardizer.uncharge_mol(m)
            self.assertEqual(esmi, Chem.MolToSmiles(nm))

    def testRemoveSGroups(self):
        mb = '''
  ACCLDraw04231812452D

 13 13  0  0  0  0  0  0  0  0999 V2000
    5.0960   -4.3327    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.4016   -4.3321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0960   -5.6718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2539   -6.3314    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2510   -2.3329    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.4056   -1.6662    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0963   -1.6662    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2510   -3.6663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5505   -6.3291    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5505   -7.6624    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.7052   -5.6625    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.4016   -5.6658    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9416   -2.3329    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  3  4  2  0  0  0  0
  1  3  1  0  0  0  0
  5  8  1  0  0  0  0
  5  7  2  0  0  0  0
  5  6  1  0  0  0  0
  1  8  2  0  0  0  0
  2  8  1  0  0  0  0
  9 12  1  0  0  0  0
  9 11  2  0  0  0  0
  9 10  1  0  0  0  0
  4 12  1  0  0  0  0
  2 12  2  0  0  0  0
  7 13  1  0  0  0  0
M  STY  2   1 DAT   2 DAT
M  SLB  2   1   1   2   2
M  SAL   1  1  10
M  SDT   1 pH
M  SDD   1     0.0000    0.0000    DR    ALL  1       6
M  SED   1 4.6
M  SAL   2  2   5   7
M  SBL   2  1   4
M  SDT   2 Stereo
M  SDD   2     0.0000    0.0000    DR    ALL  1       6
M  SED   2 E/Z unknown
M  END
'''
        m = Chem.MolFromMolBlock(mb)
        self.assertEqual(len(Chem.GetMolSubstanceGroups(m)), 2)
        nm = standardizer.remove_sgroups_from_mol(m)
        self.assertEqual(len(Chem.GetMolSubstanceGroups(nm)), 0)
        # note that the molecule is also modified in place here:
        self.assertEqual(len(Chem.GetMolSubstanceGroups(m)), 0)

    def testRemoveHs(self):
        mb = '''
  Mrv1810 04011911282D

  6  5  0  0  0  0            999 V2000
  -10.2573   -6.3835    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  -11.0823   -6.3835    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
  -10.6698   -5.6690    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.6698   -7.0979    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -11.9073   -6.3835    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  -11.0823   -7.8124    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  6  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  2  5  1  0  0  0  0
  4  6  1  0  0  0  0
M  END
'''
        m = Chem.MolFromMolBlock(mb, sanitize=False)
        self.assertEqual(m.GetNumAtoms(), 6)
        nm = standardizer.remove_hs_from_mol(m)
        self.assertEqual(nm.GetNumAtoms(), 5)
        nm = standardizer.standardize_mol(m)
        self.assertEqual(nm.GetNumAtoms(), 5)

    def testNormalize_alkoxides(self):
        mb = '''
  Mrv1810 04011911282D

  2  1  0  0  0  0            999 V2000
   -1.5000    0.0000    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
M  END
'''
        m = Chem.MolFromMolBlock(mb, sanitize=False)
        nm = standardizer.normalize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "[Na+].[OH-]")
        nm = standardizer.standardize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "[Na+].[OH-]")

        m = Chem.MolFromMolBlock(mb.replace('Na', 'K '), sanitize=False)
        nm = standardizer.normalize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "[K+].[OH-]")
        nm = standardizer.standardize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "[K+].[OH-]")

        m = Chem.MolFromMolBlock(mb.replace('Na', 'Li'), sanitize=False)
        nm = standardizer.normalize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "[Li+].[OH-]")
        nm = standardizer.standardize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "[Li+].[OH-]")

        mb = '''
  Mrv1810 04011911282D

  3  2  0  0  0  0            999 V2000
   -1.5000    0.0000    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  3  2  1  0  0  0  0
M  END
'''
        m = Chem.MolFromMolBlock(mb, sanitize=False)
        nm = standardizer.normalize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "C[O-].[Na+]")
        nm = standardizer.standardize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "C[O-].[Na+]")

        mb = '''
  Mrv1810 04171913242D

 18 18  0  0  0  0            999 V2000
   -7.0376    6.9378    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7521    6.5253    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7521    5.7003    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.0376    5.2878    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3231    5.7003    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3231    6.5253    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4665    6.9378    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.0376    7.7628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6087    6.9378    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6087    5.2878    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4665    5.2878    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.0376    4.4628    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3231    8.1753    0.0000 K   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6086    7.7628    0.0000 K   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8942    5.7003    0.0000 K   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3231    4.0503    0.0000 K   0  0  0  0  0  0  0  0  0  0  0  0
   -8.4665    4.4628    0.0000 K   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1810    6.5253    0.0000 K   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  1  6  1  0  0  0  0
  2  7  1  0  0  0  0
  1  8  1  0  0  0  0
  6  9  1  0  0  0  0
  5 10  1  0  0  0  0
  3 11  1  0  0  0  0
  4 12  1  0  0  0  0
  8 13  1  0  0  0  0
  9 14  1  0  0  0  0
 10 15  1  0  0  0  0
 12 16  1  0  0  0  0
 11 17  1  0  0  0  0
  7 18  1  0  0  0  0
M  END
'''
        m = Chem.MolFromMolBlock(mb, sanitize=False)
        nm = standardizer.normalize_mol(m)
        self.assertEqual(
            Chem.MolToSmiles(nm),
            "[K+].[K+].[K+].[K+].[K+].[K+].[O-]C1C([O-])C([O-])C([O-])C([O-])C1[O-]"
        )
        nm = standardizer.standardize_mol(m)
        self.assertEqual(
            Chem.MolToSmiles(nm),
            "[K+].[K+].[K+].[K+].[K+].[K+].[O-]C1C([O-])C([O-])C([O-])C([O-])C1[O-]"
        )

    def testNormalize_nitro(self):
        mb = '''
  Mrv1810 04031914242D

  4  3  0  0  0  0            999 V2000
   -5.5714    2.7744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8570    3.1869    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1425    2.7744    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8570    4.0119    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  4  2  0  0  0  0
  2  3  2  0  0  0  0
M  END
'''
        m = Chem.MolFromMolBlock(mb, sanitize=False)
        self.assertEqual(Chem.MolToSmiles(m), "CN(=O)=O")
        nm = standardizer.normalize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "C[N+](=O)[O-]")
        nm = standardizer.standardize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "C[N+](=O)[O-]")

    def testNormalize_sulfoxide(self):
        mb = '''
  Mrv1810 04031914252D

  4  3  0  0  0  0            999 V2000
  -10.8271    2.9774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.1126    3.3899    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3981    2.9774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.1126    4.2149    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  2  0  0  0  0
M  END
'''
        m = Chem.MolFromMolBlock(mb, sanitize=False)
        self.assertEqual(Chem.MolToSmiles(m), "CS(C)=O")
        nm = standardizer.normalize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "C[S+](C)[O-]")
        nm = standardizer.standardize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "C[S+](C)[O-]")

        mb = '''
  Mrv1810 04031914272D

  5  4  0  0  0  0            999 V2000
  -10.8271    2.9774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.1126    3.3899    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3981    2.9774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.1126    4.2149    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3981    3.8024    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  2  0  0  0  0
  2  5  2  0  0  0  0
M  END
'''
        m = Chem.MolFromMolBlock(mb, sanitize=False)
        self.assertEqual(Chem.MolToSmiles(m), "CS(C)(=O)=O")
        nm = standardizer.normalize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "CS(C)(=O)=O")
        nm = standardizer.standardize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "CS(C)(=O)=O")

        mb = '''
  Mrv1810 04031914292D

  6  5  0  0  0  0            999 V2000
  -10.8271    2.9774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.1126    3.3899    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3981    2.9774    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  -10.5637    4.0796    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.6837    3.3899    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.4995    3.9420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  2  0  0  0  0
  3  5  1  0  0  0  0
  2  6  2  0  0  0  0
M  END
'''
        m = Chem.MolFromMolBlock(mb, sanitize=False)
        self.assertEqual(Chem.MolToSmiles(m), "CNS(C)(=O)=O")
        nm = standardizer.normalize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "CNS(C)(=O)=O")
        nm = standardizer.standardize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "CNS(C)(=O)=O")

    def testNormalize1(self):
        tests = [
            ('CN(C)(C)C', 'C[N+](C)(C)C'),
            ('CN#N', 'C[N+]#N'),
            ('c1ccccc1N#N', 'c1ccccc1[N+]#N'),
            ('C=O-C', 'C=[O+]-C'),
            ('C-O-C', 'C-O-C'),
            ('O=S-C', 'O=[S+]-C'),
            ('C.Cl', 'C.Cl'),
            ('C.[Cl]', 'C.[Cl-]'),
        ]
        for smi, expected in tests:
            sp = Chem.SmilesParserParams()
            sp.removeHs = False
            em = Chem.MolFromSmiles(expected, sp)
            esmi = Chem.MolToSmiles(em)

            m = Chem.MolFromSmiles(smi, sanitize=False)
            # simulate having come from a mol file by adding coords and
            # wedging bonds:
            rdDepictor.Compute2DCoords(m)
            Chem.WedgeMolBonds(m, m.GetConformer())

            ssmi = Chem.MolToSmiles(standardizer.normalize_mol(Chem.Mol(m)))
            self.assertEqual(ssmi, esmi)

    def testUncharge_amine_hydrochlorides(self):
        mb = '''
  Mrv1810 04031914352D

  3  1  0  0  0  0            999 V2000
   -8.9098    1.3759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1953    1.7884    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   -7.1729    1.7143    0.0000 Cl  0  5  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  CHG  2   2   1   3  -1
M  END
'''
        m = Chem.MolFromMolBlock(mb, sanitize=False)
        m.UpdatePropertyCache(False)
        self.assertEqual(Chem.MolToSmiles(m), "C[NH3+].[Cl-]")
        nm = standardizer.uncharge_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "CN.Cl")
        nm = standardizer.standardize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "CN.Cl")

    def testNormalize_amide(self):
        mb = '''
  Mrv1810 04031915162D

  4  3  0  0  0  0            999 V2000
   -7.3083    2.1203    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5938    2.5328    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8793    2.1203    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5938    3.3578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  2  0  0  0  0
M  END
'''
        m = Chem.MolFromMolBlock(mb, sanitize=False)
        self.assertEqual(Chem.MolToSmiles(m), "CC(=N)O")
        nm = standardizer.normalize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "CC(N)=O")
        nm = standardizer.standardize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "CC(N)=O")

        mb = '''
  Mrv1810 04031915192D

  5  4  0  0  0  0            999 V2000
   -7.3083    2.1203    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5938    2.5328    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8793    2.1203    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.5938    3.3578    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.8793    3.7703    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  2  0  0  0  0
  4  5  1  4  0  0  0
M  END
'''
        m = Chem.MolFromMolBlock(mb, sanitize=False)
        self.assertEqual(Chem.MolToSmiles(m), "CN=C(C)O")
        nm = standardizer.normalize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "CNC(C)=O")
        nm = standardizer.standardize_mol(m)
        self.assertEqual(Chem.MolToSmiles(nm), "CNC(C)=O")

    def test_redraw_internals(self):
        m = Chem.MolFromSmiles('C1C(C1)C#CC')
        self.assertEqual(
            sorted(
                standardizer._getAtomsToOtherSide(m.GetAtomWithIdx(3),
                                                  m.GetBondBetweenAtoms(3,
                                                                        4))),
            [0, 1, 2])
        self.assertEqual(
            sorted(
                standardizer._getAtomsToOtherSide(m.GetAtomWithIdx(4),
                                                  m.GetBondBetweenAtoms(3,
                                                                        4))),
            [5])

        m = Chem.MolFromSmiles('CC#N')
        self.assertEqual(
            sorted(
                standardizer._getAtomsToOtherSide(m.GetAtomWithIdx(1),
                                                  m.GetBondBetweenAtoms(1,
                                                                        2))),
            [0])
        self.assertEqual(
            sorted(
                standardizer._getAtomsToOtherSide(m.GetAtomWithIdx(2),
                                                  m.GetBondBetweenAtoms(1,
                                                                        2))),
            [])

    def test_triple_bonds_and_allenes(self):
        ms = [x for x in Chem.SDMolSupplier(os.path.join(_test_data_dir, "odd_drawings.sdf"))]
        self.assertEqual(len(ms), 4)
        for m in ms:
            cm = standardizer.cleanup_drawing_mol(m)
            conf = cm.GetConformer()
            tmatches = m.GetSubstructMatches(Chem.MolFromSmarts('[#6]C#*'))
            amatches = m.GetSubstructMatches(Chem.MolFromSmarts('*C=C=C*'))
            self.assertTrue(len(tmatches) or len(amatches))
            for match in tmatches:
                self.assertAlmostEqual(rdMolTransforms.GetAngleRad(
                    conf, match[0], match[1], match[2]),
                                       math.pi,
                                       places=2)
            for match in amatches:
                self.assertAlmostEqual(rdMolTransforms.GetAngleRad(
                    conf, match[1], match[2], match[3]),
                                       math.pi,
                                       places=2)

    def test_cleanup_confs(self):
        mb = """
  Mrv1810 02111910063D

  2  1  0  0  0  0            999 V2000
   -0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3127    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        m = Chem.MolFromMolBlock(mb, sanitize=False, removeHs=False)
        self.assertTrue(m.GetConformer().Is3D())
        nm = standardizer.cleanup_drawing_mol(m)
        self.assertFalse(nm.GetConformer().Is3D())

        mb = """
  Mrv1810 02111910063D

  2  1  0  0  0  0            999 V2000
   -0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3127    0.4125    0.1000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        m = Chem.MolFromMolBlock(mb, sanitize=False, removeHs=False)
        self.assertTrue(m.GetConformer().Is3D())
        with self.assertRaises(ValueError):
            nm = standardizer.cleanup_drawing_mol(m)

    def test_fragment_parent1(self):
        tests = [
            ('c1cccnc1C(=O)O.[Na]', 'c1cccnc1C(=O)O'),
            ('c1cccnc1C(=O)[O-].[Na+]', 'c1cccnc1C(=O)O'),
            ('[Na].[Cl]', '[Na].[Cl]'),
            ('[Na+].[Cl-]', '[Na+].[Cl-]'),
            ('c1cccnc1[NH3+].O=C(O)C(O)C(O)C(=O)O', 'c1cccnc1N'),
            ('c1cccnc1[NH3+].O=C([O-])[C@H](O)[C@H](O)C(=O)[O-]', 'c1cccnc1N'),
            ('c1cccnc1.ClCCl', 'c1cccnc1'),
            ('c1cccnc1.ClCCl.[Na+].[Cl-].O', 'c1cccnc1'),
            ('O=C([O-])C(O)C(O)C(=O)[O-]', 'O=C(O)C(O)C(O)C(=O)O'),
            ('O=C([O-])C(O)C(O)C(=O)[O-].[Na+].[Na+]',
             'O=C([O-])C(O)C(O)C(=O)[O-].[Na+].[Na+]'),
            ('COc1cc(NCc2ccc3nc(N)nc(N)c3c2C)cc(OC)c1OC.O=CC(O)C(O)C(O)C(O)C(=O)O',
             'COc1cc(NCc2ccc3nc(N)nc(N)c3c2C)cc(OC)c1OC')
        ]
        for smi, expected in tests:
            m = Chem.MolFromSmiles(smi)
            ssmi = Chem.MolToSmiles(standardizer.get_parent_mol(m)[0])
            esmi = Chem.CanonSmiles(expected)
            self.assertEqual(ssmi, esmi)

    def test_fragment_parent2(self):
        tests = [
            ('c1cccnc1C(=O)[O-].[Na+]', 'c1cccnc1C(=O)[O-]'),
            ('[Na+].[Cl-]', '[Na+].[Cl-]'),
            ('c1cccnc1[NH3+].O=C([O-])C(O)C(O)C(=O)[O-]', 'c1cccnc1[NH3+]'),
            ('c1cccnc1[NH3+].O=C([O-])[C@H](O)[C@H](O)C(=O)[O-]',
             'c1cccnc1[NH3+]'),
            ('O=C([O-])C(O)C(O)C(=O)[O-]', 'O=C([O-])C(O)C(O)C(=O)[O-]'),
            ('O=C([O-])C(O)C(O)C(=O)[O-].[Na+].[Na+]',
             'O=C([O-])C(O)C(O)C(=O)[O-].[Na+].[Na+]'),
        ]
        for smi, expected in tests:
            m = Chem.MolFromSmiles(smi)
            ssmi = Chem.MolToSmiles(
                standardizer.get_parent_mol(m, neutralize=False)[0])
            esmi = Chem.CanonSmiles(expected)
            self.assertEqual(ssmi, esmi)
            # get_fragment_parent_mol doesn't do neutralization:
            ssmi = Chem.MolToSmiles(standardizer.get_fragment_parent_mol(m)[0])

    def test_isotopes_parent1(self):
        tests = [
            ('c1cc[13cH]nc1', 'c1cccnc1'),
            ('c1ccc([2H])nc1', 'c1cccnc1'),
            ('F[C@]([2H])(Cl)C', 'F[C@]([H])(Cl)C'),
        ]
        for smi, expected in tests:
            m = Chem.MolFromSmiles(smi)
            # simulate having come from a mol file by adding coords and
            # wedging bonds:
            rdDepictor.Compute2DCoords(m)
            Chem.WedgeMolBonds(m, m.GetConformer())
            ssmi = Chem.MolToSmiles(standardizer.get_isotope_parent_mol(m))
            sp = Chem.SmilesParserParams()
            sp.removeHs = False
            em = Chem.MolFromSmiles(expected, sp)
            esmi = Chem.MolToSmiles(em)
            self.assertEqual(ssmi, esmi)

    def test_get_parent1(self):
        tests = [
            ('c1cc[13cH]nc1.[Na]', 'c1cccnc1'),
            ('c1ccc([2H])nc1.c1ccccc1C(=O)O', 'c1cccnc1'),
            ('F[C@]([2H])(Cl)C.Cl', 'F[C@]([H])(Cl)C'),
            ('c1cccnc1.ClCCl', 'c1cccnc1'),
            ('c1cccnc1.ClCCl.[Na+].[Cl-].O', 'c1cccnc1'),
            ('O=C([O-])C(O)C(O)C(=O)[O-]', 'O=C(O)C(O)C(O)C(=O)O'),
            ('O=C([O-])C(O)C(O)C(=O)[O-].[Na+].[Na+]',
             'O=C([O-])C(O)C(O)C(=O)[O-].[Na+].[Na+]'),
            ('c1cccnc1C(=O)[O-].[Na+]', 'c1cccnc1C(=O)O'),
        ]
        for smi, expected in tests:
            m = Chem.MolFromSmiles(smi)
            # simulate having come from a mol file by adding coords and
            # wedging bonds:
            rdDepictor.Compute2DCoords(m)
            Chem.WedgeMolBonds(m, m.GetConformer())
            ssmi = Chem.MolToSmiles(standardizer.get_parent_mol(m)[0])
            sp = Chem.SmilesParserParams()
            sp.removeHs = False
            em = Chem.MolFromSmiles(expected, sp)
            esmi = Chem.MolToSmiles(em)
            print(smi)
            self.assertEqual(ssmi, esmi)

            # make sure we can also work from mol blocks:
            imb = Chem.MolToMolBlock(m)
            smb, exclude = standardizer.get_parent_molblock(imb)
            sm = Chem.MolFromMolBlock(smb, removeHs=False)
            ssmi = Chem.MolToSmiles(sm)
            self.assertEqual(ssmi, esmi)

    def test_exclude(self):
        molb = '''
  Mrv1810 04251917202D          

  3  2  0  0  0  0            999 V2000
   -0.8797    1.4887    0.0000 Mo  0  0  0  0  0  0  0  0  0  0  0  0
   -0.9229    0.6649    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -1.7002    1.5750    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
'''
        mol = Chem.MolFromMolBlock(molb, sanitize=False, removeHs=False)
        omol = standardizer.standardize_mol(mol)
        # should have returned the same object
        self.assertTrue(mol is omol)

        molb = '''
  Mrv1810 04251917202D          

  3  2  0  0  0  0            999 V2000
   -0.8797    1.4887    0.0000 Fe  0  0  0  0  0  0  0  0  0  0  0  0
   -0.9229    0.6649    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -1.7002    1.5750    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
'''
        mol = Chem.MolFromMolBlock(molb, sanitize=False, removeHs=False)
        omol = standardizer.standardize_mol(mol)
        # We no longer standardize this
        self.assertTrue(mol is omol)

    def test_keep_hs(self):
        sdf = '''
  SciTegic05151911262D

 25 28  0  0  0  0            999 V2000
    2.3968   -0.6245    0.0000 C   0  0  2  0  0  0
    3.1901   -0.6695    0.0000 C   0  0
    2.0086    0.0788    0.0000 C   0  0
    3.6121   -1.3784    0.0000 C   0  0
    1.9692   -1.3109    0.0000 C   0  0
    2.7231   -1.9017    0.0000 O   0  0
    2.4193    0.7764    0.0000 C   0  0
    3.6459    0.0113    0.0000 C   0  0
    1.1590    0.1238    0.0000 C   0  0
    3.2520    0.7483    0.0000 C   0  0
    2.0424    1.4853    0.0000 N   0  0
    2.8807    0.0000    0.0000 C   0  0
    0.7145   -0.5851    0.0000 C   0  0
    1.1253   -1.2884    0.0000 C   0  0  2  0  0  0
    4.4617   -1.3897    0.0000 C   0  0
    2.9088    1.4853    0.0000 C   0  0
    4.4842    0.0113    0.0000 C   0  0
    6.3127   -0.1041    0.0000 Cl  0  0
    4.9062   -0.6977    0.0000 C   0  0
    4.7544   -1.9981    0.0000 O   0  0
    0.7497   -1.8495    0.0000 O   0  0
    1.4729    1.8481    0.0000 C   0  0
    1.7424    0.5744    0.0000 H   0  0
    2.7172    1.2538    0.0000 H   0  0
    1.7786   -1.8403    0.0000 H   0  0
  1  2  1  0
  1  3  1  0
  2  4  1  0
  1  5  1  0
  4  6  1  0
  5  6  1  0
  3  7  1  0
  2  8  2  0
  3  9  1  0
  7 10  1  0
  8 10  1  0
  7 11  1  0
  1 12  1  6
  9 13  2  0
  5 14  1  0
 13 14  1  0
  4 15  2  0
 11 16  1  0
 12 16  1  0
  8 17  1  0
 15 19  1  0
 17 19  2  0
 15 20  1  0
 14 21  1  1
 11 22  1  0
  3 23  1  0
  7 24  1  0
  5 25  1  0
M  END
> <mol_id>
1

$$$$

  SciTegic05151911262D

 25 28  0  0  0  0            999 V2000
    4.2600   -1.1100    0.0000 C   0  0  2  0  0  0
    5.6700   -1.1900    0.0000 C   0  0
    3.5700    0.1400    0.0000 C   0  0  2  0  0  0
    6.4200   -2.4500    0.0000 C   0  0
    3.5000   -2.3300    0.0000 C   0  0  1  0  0  0
    4.8400   -3.3800    0.0000 O   0  0
    4.3000    1.3800    0.0000 C   0  0  1  0  0  0
    6.4800    0.0200    0.0000 C   0  0
    2.0600    0.2200    0.0000 C   0  0
    5.7800    1.3300    0.0000 C   0  0
    3.6300    2.6400    0.0000 N   0  0
    5.1200    0.0000    0.0000 C   0  0
    1.2700   -1.0400    0.0000 C   0  0
    2.0000   -2.2900    0.0000 C   0  0  2  0  0  0
    7.9300   -2.4700    0.0000 C   0  0
    5.1700    2.6400    0.0000 C   0  0
    7.9700    0.0200    0.0000 C   0  0
   11.2200   -0.1850    0.0000 Cl  0  0
    8.7200   -1.2400    0.0000 C   0  0
    8.4503   -3.5513    0.0000 O   0  0
    1.3325   -3.2873    0.0000 O   0  0
    2.6179    3.2848    0.0000 C   0  0
    3.0968    1.0210    0.0000 H   0  0
    4.8294    2.2284    0.0000 H   0  0
    3.1612   -3.2709    0.0000 H   0  0
  1  2  1  0
  1  3  1  0
  2  4  1  0
  1  5  1  0
  4  6  1  0
  5  6  1  0
  3  7  1  0
  2  8  2  0
  3  9  1  0
  7 10  1  0
  8 10  1  0
  7 11  1  0
  1 12  1  6
  9 13  2  0
  5 14  1  0
 13 14  1  0
  4 15  2  0
 11 16  1  0
 12 16  1  0
  8 17  1  0
 15 19  1  0
 17 19  2  0
 15 20  1  0
 14 21  1  1
 11 22  1  0
  3 23  1  6
  7 24  1  1
  5 25  1  6
M  END
> <mol_id>
2

$$$$

  SciTegic05151911262D

 32 35  0  0  0  0            999 V2000
   -2.3925    1.1451    0.0000 C   0  0
   -3.1994    1.3166    0.0000 C   0  0
   -3.7515    0.7035    0.0000 C   0  0  1  0  0  0
   -3.4965   -0.0811    0.0000 C   0  0
   -2.6896   -0.2526    0.0000 C   0  0
   -2.4346   -1.0373    0.0000 C   0  0
   -1.6276   -1.2088    0.0000 C   0  0
   -1.0756   -0.5957    0.0000 C   0  0  2  0  0  0
   -1.3306    0.1889    0.0000 C   0  0  1  0  0  0
   -2.1375    0.3605    0.0000 C   0  0  1  0  0  0
   -0.7785    0.8020    0.0000 C   0  0
    0.0285    0.6305    0.0000 C   0  0
    0.2834   -0.1541    0.0000 C   0  0  1  0  0  0
   -0.2686   -0.7672    0.0000 C   0  0  1  0  0  0
    0.1439   -1.4817    0.0000 C   0  0
    0.9508   -1.3102    0.0000 C   0  0
    1.0371   -0.4897    0.0000 C   0  0  2  0  0  0
    0.8965    0.3979    0.0000 C   0  0
   -1.5855    0.9736    0.0000 C   0  0
    1.7515   -0.0772    0.0000 C   0  0  2  0  0  0
    1.7515    0.7478    0.0000 C   0  0
    2.4660   -0.4897    0.0000 C   0  0
    3.1805   -0.0772    0.0000 C   0  0
    3.8950   -0.4897    0.0000 C   0  0
    4.6094   -0.0772    0.0000 C   0  0
    5.3239   -0.4897    0.0000 C   0  0
    4.6094    0.7478    0.0000 C   0  0
   -4.5584    0.8750    0.0000 O   0  0
    1.8758   -1.0344    0.0000 H   0  0
   -0.8564   -1.5762    0.0000 H   0  0
   -1.6395    1.1400    0.0000 H   0  0
   -0.7665   -1.5467    0.0000 H   0  0
  1  2  1  0
  1 10  1  0
  2  3  1  0
  3  4  1  0
  3 28  1  1
  4  5  1  0
  5  6  2  0
 10  5  1  0
  6  7  1  0
  8  7  1  0
  8  9  1  0
 14  8  1  0
 10  9  1  0
  9 11  1  0
 10 19  1  1
 11 12  1  0
 12 13  1  0
 13 14  1  0
 13 17  1  0
 13 18  1  1
 14 15  1  0
 15 16  1  0
 16 17  1  0
 17 20  1  0
 20 21  1  6
 20 22  1  0
 22 23  1  0
 23 24  1  0
 24 25  1  0
 25 26  1  0
 25 27  1  0
 17 29  1  6
 14 30  1  6
  9 31  1  6
  8 32  1  1
M  END
> <mol_id>
3

$$$$

  SciTegic05151911262D

  5  4  0  0  0  0            999 V2000
    7.4839   -1.1196    0.0000 P   0  0
    7.4839   -0.2946    0.0000 O   0  0
    6.6589   -1.1196    0.0000 O   0  0
    8.0673   -1.7030    0.0000 H   0  0
    7.2704   -1.9165    0.0000 H   0  0
  1  2  2  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  END
> <mol_id>
8

$$$$
1266
  -OEChem-01301907122D

 21 22  0     1  0  0  0  0  0999 V2000
    7.2800    4.2017    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   12.0627    1.4154    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1479    3.0026    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.8373    2.5483    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    9.8506    0.2719    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.4268    4.2160    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    8.5018    2.1058    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5101    0.6989    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2906    2.7986    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.6711    1.4150    0.0000 C  -1  0  0  0  0  0  0  0  0  0  0  0
    7.2937    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0854    2.0997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8534    4.2088    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
    6.0835    0.6989    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0668    4.9001    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6400    4.9072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2134    4.9143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    4.2229    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5079    6.1275    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9046    6.1275    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5446    3.0026    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  9  1  0  0  0  0
  1 15  1  0  0  0  0
  2 10  2  0  0  0  0
 13  3  1  1  0  0  0
  4  7  1  0  0  0  0
  4 10  1  0  0  0  0
  5  8  1  0  0  0  0
  5 10  1  0  0  0  0
  6 16  1  0  0  0  0
  6 17  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  2  0  0  0  0
  8 11  2  0  0  0  0
  9 12  1  0  0  0  0
 11 14  1  0  0  0  0
 12 14  2  0  0  0  0
 13 15  1  0  0  0  0
 13 16  1  6  0  0  0
 13 21  1  0  0  0  0
 17 18  1  0  0  0  0
 17 19  1  0  0  0  0
 17 20  1  0  0  0  0
M  ISO  1  10  11
M  END
        '''
        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf, removeHs=False, sanitize=False)
        ms = [x for x in suppl]
        self.assertEqual(len(ms), 5)
        q = rdqueries.AtomNumEqualsQueryAtom(1)
        for m in ms:
            self.assertTrue(len(m.GetAtomsMatchingQuery(q)) >= 1)
            nm = standardizer.remove_hs_from_mol(m)
            self.assertEqual(m.GetNumAtoms(), nm.GetNumAtoms())

        # and examples where Hs should be removed:
        sdf = '''7625701
  -OEChem-01301907452D

 40 42  0     0  0  0  0  0  0999 V2000
   -1.7436    6.7714    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -3.1026    2.9970    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4510    5.2394    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4394    1.8840    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000    2.0104    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   -0.8675    1.5027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8675    1.5027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8675    0.4975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8675    0.4975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1275    3.3488    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7718    4.1135    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7641    3.9423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4350    5.0613    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4115    4.7135    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0827    5.8306    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0673    5.6524    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3978    4.5394    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7443    3.5947    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0967    2.8235    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0397    5.3062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4832    2.5840    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6058    1.3710    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1268    2.2065    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1268    2.2065    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6058    1.3710    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1227   -0.2078    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6062    0.6270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6062    0.6270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1227   -0.2078    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4815   -0.5750    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4815   -0.5750    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7011    2.8655    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5539    3.8320    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5500    6.2264    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4830    3.4652    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1978    5.9454    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6148    4.8248    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4646    5.7876    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5211    5.8813    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1 16  1  0  0  0  0
  2 13  1  0  0  0  0
  2 20  1  0  0  0  0
  3 14  1  0  0  0  0
  3 37  1  0  0  0  0
  4 20  2  0  0  0  0
  5  6  1  0  0  0  0
  5  7  1  0  0  0  0
  5 11  1  0  0  0  0
  5 22  1  0  0  0  0
  6  8  1  0  0  0  0
  6 23  1  0  0  0  0
  6 24  1  0  0  0  0
  7  9  1  0  0  0  0
  7 25  1  0  0  0  0
  7 26  1  0  0  0  0
  8 10  1  0  0  0  0
  8 27  1  0  0  0  0
  8 28  1  0  0  0  0
  9 10  1  0  0  0  0
  9 29  1  0  0  0  0
  9 30  1  0  0  0  0
 10 31  1  0  0  0  0
 10 32  1  0  0  0  0
 11 12  1  0  0  0  0
 11 33  1  0  0  0  0
 11 34  1  0  0  0  0
 12 13  1  0  0  0  0
 12 14  2  0  0  0  0
 13 15  2  0  0  0  0
 14 16  1  0  0  0  0
 15 17  1  0  0  0  0
 15 18  1  0  0  0  0
 16 17  2  0  0  0  0
 17 35  1  0  0  0  0
 18 19  2  0  0  0  0
 18 21  1  0  0  0  0
 19 20  1  0  0  0  0
 19 36  1  0  0  0  0
 21 38  1  0  0  0  0
 21 39  1  0  0  0  0
 21 40  1  0  0  0  0
M  CHG  1   5   1
M  END
$$$$

  Mrv1810 05161910292D          

  5  4  0  0  0  0            999 V2000
   -2.3910    2.3684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5660    2.3684    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   -1.1535    1.6540    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1535    3.0829    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9785    3.0829    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  2  5  1  0  0  0  0
M  CHG  1   2   1
M  END
$$$$

  Mrv1810 05161910302D          

  5  4  0  0  0  0            999 V2000
   -2.3910    2.3684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5660    2.3684    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   -1.1535    1.6540    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1535    3.0829    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9785    3.0829    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  2  5  1  0  0  0  0
M  CHG  1   2   1
M  END
$$$
'''
        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf, removeHs=False, sanitize=False)
        ms = [x for x in suppl]
        self.assertEqual(len(ms), 3)
        q = rdqueries.AtomNumEqualsQueryAtom(1)
        for m in ms:
            self.assertTrue(len(m.GetAtomsMatchingQuery(q)) >= 1)
            nm = standardizer.remove_hs_from_mol(m)
            self.assertTrue(m.GetNumAtoms() > nm.GetNumAtoms())
            self.assertEqual(len(nm.GetAtomsMatchingQuery(q)), 0)

    def test_charge_h_interaction(self):
        sdf = '''
  Mrv1810 05161910292D

  5  4  0  0  0  0            999 V2000
   -2.3910    2.3684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5660    2.3684    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   -1.1535    1.6540    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1535    3.0829    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9785    3.0829    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  2  5  1  0  0  0  0
M  CHG  1   2   1
M  END
$$$$

  Mrv1810 05161910302D

  5  4  0  0  0  0            999 V2000
   -2.3910    2.3684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5660    2.3684    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   -1.1535    1.6540    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1535    3.0829    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9785    3.0829    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  2  5  1  0  0  0  0
M  CHG  1   2   1
M  END
$$$$

  Mrv1810 05161911452D          

 15 17  0  0  0  0            999 V2000
    3.1497    3.8099    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5951    1.6260    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    3.0372    1.1315    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    2.2137    1.0800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5317    2.5737    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0895    3.0682    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0903    2.1620    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5832    1.7502    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5436    2.4495    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8923    1.9499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9130    3.1197    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4291    0.8251    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6987    0.4126    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1035    1.3855    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.2778    0.6399    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1 11  2  0  0  0  0
  2  4  1  0  0  0  0
  2  7  1  0  0  0  0
  2  9  1  0  0  0  0
  2 14  1  0  0  0  0
  3  4  1  0  0  0  0
  3  8  1  0  0  0  0
  3 10  1  0  0  0  0
  3 15  1  0  0  0  0
  4 12  1  0  0  0  0
  4 13  1  0  0  0  0
  5  7  1  0  0  0  0
  5  8  1  0  0  0  0
  5 11  1  0  0  0  0
  6  9  1  0  0  0  0
  6 10  1  0  0  0  0
  6 11  1  0  0  0  0
M  CHG  2   2   1   3   1
M  END
$$$$
'''
        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf, removeHs=False, sanitize=False)
        ms = [x for x in suppl]
        self.assertEqual(len(ms), 3)
        q = rdqueries.AtomNumEqualsQueryAtom(1)
        for m in ms:
            self.assertTrue(len(m.GetAtomsMatchingQuery(q)) >= 1)
            print('--------------------------')
            nm = standardizer.standardize_mol(m)
            print('--------------------------')
            self.assertTrue(m.GetNumAtoms() > nm.GetNumAtoms())
            self.assertEqual(len(nm.GetAtomsMatchingQuery(q)), 0)
            self.assertEqual(Chem.GetFormalCharge(nm), 0)

    def test_preserve_wedging(self):
        molb = '''4116
  -OEChem-01301907122D

 27 31  0     1  0  0  0  0  0999 V2000
   20.9488  -13.3575    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   20.4708  -15.1657    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   26.6044  -15.8799    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   17.6777  -16.5018    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   23.2715  -19.1019    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   20.9604  -19.0322    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   20.9488  -15.8799    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   20.9488  -14.6130    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   23.1431  -14.6130    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   23.1431  -15.8799    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   22.0489  -16.5134    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   19.8547  -16.5192    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   19.4227  -15.2119    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   24.2488  -16.5134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   19.8547  -13.9851    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   25.3546  -15.8857    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   23.1374  -17.5155    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   24.2431  -13.9679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   18.7662  -15.8799    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   22.0489  -18.4103    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   25.3488  -17.5155    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   18.7662  -14.6130    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   25.3488  -14.6014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   19.6389  -17.7690    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   26.2358  -18.4025    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   20.9488  -17.1997    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   23.1374  -13.3516    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  8  1  1  6  0  0  0
  1 13  1  0  0  0  0
  2 13  2  0  0  0  0
 16  3  1  6  0  0  0
 19  4  1  1  0  0  0
  5 20  2  0  0  0  0
  6 20  1  0  0  0  0
  7  8  1  0  0  0  0
  7 11  1  0  0  0  0
  7 12  1  0  0  0  0
  7 26  1  1  0  0  0
  8  9  1  0  0  0  0
  8 15  1  0  0  0  0
  9 10  1  0  0  0  0
  9 18  1  0  0  0  0
  9 27  1  1  0  0  0
 10 11  1  0  0  0  0
 10 14  1  0  0  0  0
 10 17  1  1  0  0  0
 11 20  1  1  0  0  0
 12 13  1  6  0  0  0
 12 19  1  0  0  0  0
 12 24  1  1  0  0  0
 14 16  1  0  0  0  0
 15 22  1  0  0  0  0
 16 21  1  1  0  0  0
 16 23  1  0  0  0  0
 17 21  1  0  0  0  0
 18 23  1  0  0  0  0
 19 22  1  0  0  0  0
 21 25  2  0  0  0  0
M  END'''
        bonds = '''  8  1  1  6  0  0  0
  1 13  1  0  0  0  0
  2 13  2  0  0  0  0
 16  3  1  6  0  0  0
 19  4  1  1  0  0  0
  5 20  2  0  0  0  0
  6 20  1  0  0  0  0
  7  8  1  0  0  0  0
  7 11  1  0  0  0  0
  7 12  1  0  0  0  0
  7 26  1  1  0  0  0
  8  9  1  0  0  0  0
  8 15  1  0  0  0  0
  9 10  1  0  0  0  0
  9 18  1  0  0  0  0
  9 27  1  1  0  0  0
 10 11  1  0  0  0  0
 10 14  1  0  0  0  0
 10 17  1  1  0  0  0
 11 20  1  1  0  0  0
 12 13  1  6  0  0  0
 12 19  1  0  0  0  0
 12 24  1  1  0  0  0
 14 16  1  0  0  0  0
 15 22  1  0  0  0  0
 16 21  1  1  0  0  0
 16 23  1  0  0  0  0
 17 21  1  0  0  0  0
 18 23  1  0  0  0  0
 19 22  1  0  0  0  0
 21 25  2  0  0  0  0
'''.replace('  0  0  0', '')
        omolb = standardizer.standardize_molblock(molb)
        # should have returned the same object
        self.assertTrue(omolb.find(bonds) > 0)

    def test_oxaniacic_acid(self):
        molb = '''
  Mrv1810 05161913202D          

 10 10  0  0  1  0            999 V2000
    0.0000    2.4835    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -2.1445    0.4084    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4283   -0.8281    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.6585    0.0000 N   0  2  0  0  0  0  0  0  0  0  0  0
   -0.7157    0.4104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7157    1.2397    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7157    0.4104    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7157    1.2397    0.0000 C   0  5  0  0  0  0  0  0  0  0  0  0
   -1.4295   -0.0031    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  4  1  0  0  0  0
  2 10  1  0  0  0  0
  3 10  2  0  0  0  0
  4  6  1  0  0  0  0
  5  6  2  0  0  0  0
  5  7  1  0  0  0  0
  5 10  1  0  0  0  0
  7  8  2  0  0  0  0
  8  9  1  0  0  0  0
  4  9  1  0  0  0  0
M  CHG  3   1  -1   4   2   9  -1
M  END
'''
        omolb = standardizer.standardize_molblock(molb)
        omol = Chem.MolFromMolBlock(omolb, sanitize=False)
        self.assertEqual(omol.GetAtomWithIdx(1).GetAtomicNum(), 7)
        self.assertEqual(omol.GetAtomWithIdx(1).GetFormalCharge(), 1)
        self.assertEqual(omol.GetAtomWithIdx(0).GetAtomicNum(), 6)
        self.assertEqual(omol.GetAtomWithIdx(0).GetFormalCharge(), 0)
        Chem.SanitizeMol(omol)
        self.assertEqual(Chem.MolToSmiles(omol), 'O=C(O)c1ccc[n+]([O-])c1')
        molb = '''7639098
  -OEChem-01301907482D

 31 33  0     1  0  0  0  0  0999 V2000
    5.1292    2.1857    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -3.9447   -3.9438    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.4993    2.5426    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    1.3133    0.9518    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.5008    1.5426    0.0000 N   0  2  0  0  0  0  0  0  0  0  0  0
    1.0015    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5888   -0.8082    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5832   -0.7024    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.2648    1.2595    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4712    2.2380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0051    0.5871    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3065    0.9518    0.0000 C   0  5  3  0  0  0  0  0  0  0  0  0
    3.4277    2.5473    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9615    0.8964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1777    1.8781    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1721   -1.5107    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7639   -2.4236    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1662   -1.4019    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3558   -3.2360    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7581   -2.2144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3559   -3.1356    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9139    2.7398    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.8470   -0.1460    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0199    1.1833    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5825    3.2811    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.5167    0.3922    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0180   -2.5019    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4693   -0.7159    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0496   -3.9207    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5036   -2.1328    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1 16  1  0  0  0  0
  2 22  1  0  0  0  0
  3  5  1  0  0  0  0
  4  5  1  0  0  0  0
  4  6  1  0  0  0  0
  4  9  1  0  0  0  0
  5 13  1  0  0  0  0
  6 12  2  0  0  0  0
  7  8  2  3  0  0  0
  7 12  1  0  0  0  0
  8 17  1  0  0  0  0
  9 10  2  0  0  0  0
  9 11  1  0  0  0  0
 10 14  1  0  0  0  0
 10 23  1  0  0  0  0
 11 15  2  0  0  0  0
 11 24  1  0  0  0  0
 12 13  1  0  0  0  0
 13 25  1  0  0  0  0
 14 16  2  0  0  0  0
 14 26  1  0  0  0  0
 15 16  1  0  0  0  0
 15 27  1  0  0  0  0
 17 18  2  0  0  0  0
 17 19  1  0  0  0  0
 18 20  1  0  0  0  0
 18 28  1  0  0  0  0
 19 21  2  0  0  0  0
 19 29  1  0  0  0  0
 20 22  2  0  0  0  0
 20 30  1  0  0  0  0
 21 22  1  0  0  0  0
 21 31  1  0  0  0  0
M  CHG  3   3  -1   5   2  13  -1
M  END
'''
        omolb = standardizer.standardize_molblock(molb)
        omol = Chem.MolFromMolBlock(omolb, sanitize=False)
        self.assertEqual(omol.GetAtomWithIdx(1).GetAtomicNum(), 7)
        self.assertEqual(omol.GetAtomWithIdx(1).GetFormalCharge(), 1)
        self.assertEqual(omol.GetAtomWithIdx(0).GetAtomicNum(), 6)
        self.assertEqual(omol.GetAtomWithIdx(0).GetFormalCharge(), 0)
        Chem.SanitizeMol(omol)
        self.assertEqual(Chem.MolToSmiles(omol),
                         '[O-][n+]1cc(N=Nc2ccc(Cl)cc2)nn1-c1ccc(Cl)cc1')

        molb = '''
  Mrv1810 05161914502D          

 10 10  0  0  1  0            999 V2000
    0.3874   -1.5568    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -1.7571   -3.6319    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0408   -4.8684    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3874   -2.3818    0.0000 N   0  2  0  0  0  0  0  0  0  0  0  0
   -0.3282   -3.6299    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3282   -2.8006    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3874   -4.0403    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1031   -3.6299    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1031   -2.8006    0.0000 N   0  5  0  0  0  0  0  0  0  0  0  0
   -1.0421   -4.0434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  4  1  0  0  0  0
  2 10  1  0  0  0  0
  3 10  2  0  0  0  0
  4  6  1  0  0  0  0
  5  6  2  0  0  0  0
  5  7  1  0  0  0  0
  5 10  1  0  0  0  0
  7  8  2  0  0  0  0
  8  9  1  0  0  0  0
  4  9  1  0  0  0  0
M  CHG  3   1  -1   4   2   9  -1
M  END
'''
        omolb = standardizer.standardize_molblock(molb)
        print(omolb)
        omol = Chem.MolFromMolBlock(omolb, sanitize=False)
        self.assertEqual(omol.GetAtomWithIdx(1).GetAtomicNum(), 7)
        self.assertEqual(omol.GetAtomWithIdx(1).GetFormalCharge(), 1)
        self.assertEqual(omol.GetAtomWithIdx(0).GetAtomicNum(), 7)
        self.assertEqual(omol.GetAtomWithIdx(0).GetFormalCharge(), 0)
        Chem.SanitizeMol(omol)
        self.assertEqual(Chem.MolToSmiles(omol), 'O=C(O)c1ccn[n+]([O-])c1')

    def testAmidePatternProblem(self):
        molb = '''
  Mrv0541 07061712112D          

  8  8  0  0  0  0            999 V2000
   12.1729   -9.6736    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.7643   -8.9645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.9471   -8.9645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.5385   -9.6736    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.9471  -10.3827    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.7643  -10.3827    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   12.1729  -11.0876    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   12.9901   -9.6736    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  1  6  2  0  0  0  0
  6  7  1  0  0  0  0
  1  8  1  0  0  0  0
M  CHG  2   6   1   7  -1
M  END
'''
        omolb = standardizer.standardize_molblock(molb)
        nm = Chem.MolFromMolBlock(omolb)
        smi = Chem.MolToSmiles(nm)
        self.assertEqual(smi, "[O-][n+]1ccccc1O")

    def testSulfoxideProblem(self):
        molb = '''
  Mrv1810 08151906432D          

  9  9  0  0  0  0            999 V2000
   -8.3910    0.1483    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1054   -0.2642    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1054   -1.0892    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3910   -1.5017    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6765   -1.0892    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6765   -0.2642    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9620    0.1483    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9620    0.9733    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2476   -0.2642    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  1  6  2  0  0  0  0
  6  7  1  0  0  0  0
  7  8  2  0  0  0  0
  7  9  1  0  0  0  0
M  END
'''
        omolb = standardizer.standardize_molblock(molb)
        nm = Chem.MolFromMolBlock(omolb)
        smi = Chem.MolToSmiles(nm)
        self.assertEqual(smi, "C[S+]([O-])c1ccccc1")
        molb = '''
  Mrv1810 08151906432D          

  9  9  0  0  0  0            999 V2000
   -8.3910    0.1483    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1054   -0.2642    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1054   -1.0892    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3910   -1.5017    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6765   -1.0892    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.6765   -0.2642    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9620    0.1483    0.0000 S   0  3  0  0  0  0  0  0  0  0  0  0
   -6.9620    0.9733    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2476   -0.2642    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  1  6  2  0  0  0  0
  6  7  1  0  0  0  0
  7  8  2  0  0  0  0
  7  9  1  0  0  0  0
M  CHG  1   7   1
M  END
'''
        omolb = standardizer.standardize_molblock(molb)
        nm = Chem.MolFromMolBlock(omolb)
        smi = Chem.MolToSmiles(nm)
        self.assertEqual(smi, "C[S+]([O-])c1ccccc1")

    def testValencePropagation(self):
        ''' another example found in testing '''
        # first case:
        molb = '''atomic C


  1  0  0  0  0  0            999 V2000
   -0.3958   -0.0542    0.0000 C   0  0  0  0  0 15
M  END'''
        omolb = standardizer.standardize_molblock(molb)
        # print(omolb)
        self.assertNotEqual(omolb.find('0  0 15'), -1)
        self.assertEqual(Chem.MolBlockToInchi(omolb), 'InChI=1S/C')

        # second case:
        molb = '''H3PO2


  3  2  0  0  0  0            999 V2000
    0.2667   -0.4167    0.0000 P   0  0  0  0  0  5
    0.2667    1.1083    0.0000 O   0  0
   -1.0958   -1.0042    0.0000 O   0  0
  2  1  2  0
  3  1  1  0
M  END'''
        omolb = standardizer.standardize_molblock(molb)
        # print(omolb)
        self.assertNotEqual(omolb.find('0  0  5'), -1)
        self.assertEqual(Chem.MolBlockToInchi(omolb),
                         'InChI=1S/H3O2P/c1-3-2/h3H2,(H,1,2)')

    def testDoubleBondStereoProblem1(self):
        ''' a problem found in testing '''
        # simplified version
        molb = '''
  Mrv1810 09181912182D          

 10  9  0  0  0  0            999 V2000
   -3.8346    1.6241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1201    2.0366    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4056    1.6241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5491    2.0366    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2635    1.6241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5491    2.8616    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2635    2.4491    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6912    2.0366    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4056    0.7991    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6912    1.2116    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  4  1  1  1  0  0  0
  4  5  1  0  0  0  0
  4  6  1  0  0  0  0
  4  7  1  0  0  0  0
  3  8  1  0  0  0  0
  3  9  1  0  0  0  0
  3 10  1  0  0  0  0
  3  2  1  1  0  0  0
M  END
'''
        omolb = standardizer.standardize_molblock(molb)
        # print(omolb)
        self.assertEqual(omolb.find('1  2  2  3'), -1)
        self.assertNotEqual(omolb.find('4  5  1  6'), -1)
        self.assertNotEqual(omolb.find('3  9  1  6'), -1)
        self.assertEqual(
            Chem.MolBlockToInchi(omolb),
            'InChI=1S/C6H10F2O2/c1-5(7,9)3-4-6(2,8)10/h3-4,9-10H,1-2H3/b4-3+/t5-,6+'
        )

        # the actual thing found in testing:
        molb = '''
          11280715002D 1   1.00000     0.00000     0

 29 32  0     1  0            999 V2000
    5.0000   -1.5375    0.0000 C   0  0  2  0  0  0           0  0  0
    4.2167   -1.2792    0.0000 C   0  0  0  0  0  0           0  0  0
    5.7167   -2.7750    0.0000 C   0  0  1  0  0  0           0  0  0
    5.0000   -2.3625    0.0000 C   0  0  1  0  0  0           0  0  0
    3.7292   -1.9542    0.0000 O   0  0  0  0  0  0           0  0  0
    4.2167   -2.6167    0.0000 C   0  0  1  0  0  0           0  0  0
    5.7167   -3.6000    0.0000 C   0  0  0  0  0  0           0  0  0
    5.7125   -5.2417    0.0000 N   0  0  3  0  0  0           0  0  0
    5.7167   -1.1250    0.0000 C   0  0  0  0  0  0           0  0  0
    6.4292   -2.3625    0.0000 C   0  0  1  0  0  0           0  0  0
    6.4292   -1.5375    0.0000 C   0  0  2  0  0  0           0  0  0
    5.0000   -4.0042    0.0000 C   0  0  0  0  0  0           0  0  0
    4.9917   -4.8292    0.0000 C   0  0  1  0  0  0           0  0  0
    3.9625   -0.4917    0.0000 O   0  0  0  0  0  0           0  0  0
    5.7042   -6.0792    0.0000 C   0  0  2  0  0  0           0  0  0
    6.4292   -4.8250    0.0000 C   0  0  0  0  0  0           0  0  0
    3.9625   -3.4042    0.0000 C   0  0  0  0  0  0           0  0  0
    7.1417   -2.7750    0.0000 C   0  0  0  0  0  0           0  0  0
    7.1417   -1.1250    0.0000 C   0  0  0  0  0  0           0  0  0
    4.2792   -6.0792    0.0000 C   0  0  0  0  0  0           0  0  0
    4.2792   -5.2500    0.0000 C   0  0  0  0  0  0           0  0  0
    4.9917   -6.4875    0.0000 C   0  0  0  0  0  0           0  0  0
    6.4292   -6.4917    0.0000 C   0  0  0  0  0  0           0  0  0
    7.8542   -2.3625    0.0000 C   0  0  0  0  0  0           0  0  0
    7.8542   -1.5375    0.0000 C   0  0  0  0  0  0           0  0  0
    5.0000   -0.7125    0.0000 H   0  0  0  0  0  0           0  0  0
    6.4292   -3.1875    0.0000 H   0  0  0  0  0  0           0  0  0
    6.4292   -0.7042    0.0000 H   0  0  0  0  0  0           0  0  0
    5.0000   -3.1792    0.0000 H   0  0  0  0  0  0           0  0  0
  2  1  1  0     0  0
  3  4  1  0     0  0
  4  1  1  0     0  0
  5  2  1  0     0  0
  6  4  1  0     0  0
  3  7  1  6     0  0
  8 13  1  0     0  0
  9  1  1  0     0  0
 10 11  1  0     0  0
 11  9  1  0     0  0
 12  7  2  0     0  0
 13 12  1  6     0  0
 14  2  2  0     0  0
 15  8  1  0     0  0
 16  8  1  0     0  0
  6 17  1  6     0  0
 18 10  1  0     0  0
 19 11  1  0     0  0
 20 21  1  0     0  0
 21 13  1  0     0  0
 22 20  1  0     0  0
 15 23  1  1     0  0
 24 25  1  0     0  0
 25 19  1  0     0  0
  1 26  1  6     0  0
 10 27  1  6     0  0
 11 28  1  6     0  0
  4 29  1  6     0  0
  5  6  1  0     0  0
  3 10  1  0     0  0
 18 24  1  0     0  0
 22 15  1  0     0  0
M  END'''
        omolb = standardizer.standardize_molblock(molb)
        # print(omolb)
        self.assertEqual(omolb.find('12  7  2  3'), -1)
        self.assertNotEqual(omolb.find('3  7  1  6'), -1)
        self.assertNotEqual(omolb.find('13 12  1  6'), -1)
        self.assertEqual(
            Chem.MolBlockToInchi(omolb),
            'InChI=1S/C22H35NO2/c1-14-7-6-9-17(23(14)3)11-12-19-18-10-5-4-8-16(18)13-20-21(19)15(2)25-22(20)24/h11-12,14-21H,4-10,13H2,1-3H3/b12-11+/t14-,15-,16+,17+,18+,19-,20-,21+/m0/s1'
        )
        molb = '''
  Mrv1810 09191906362D          

  5  4  0  0  0  0            999 V2000
  -13.1335    2.8866    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
  -12.3085    2.8866    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.9584    2.8866    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.8961    2.1762    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -14.3708    2.1762    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  2  0  0  0  0
  3  1  2  0  0  0  0
  4  2  1  0  0  0  0
  3  5  1  1  0  0  0
M  END
'''
        omolb = standardizer.standardize_molblock(molb)
        # print(omolb)
        self.assertNotEqual(omolb.find('2  1  2  0'), -1)
        self.assertNotEqual(omolb.find('3  1  2  0'), -1)
        self.assertNotEqual(omolb.find('3  5  1  1'), -1)
        self.assertEqual(Chem.MolBlockToInchi(omolb),
                         'InChI=1S/C5H8/c1-3-5-4-2/h3-4H,1-2H3/t5-/m1/s1')

        # another example found in testing:
        molb = '''
  SciTegic04201611242D

 49 52  0  0  0  0            999 V2000
    2.6000    0.0000    0.0000 N   0  3
    1.3000    0.7500    0.0000 C   0  0
    3.9000    0.7500    0.0000 C   0  0
    2.6031   -1.5008    0.0000 C   0  0
    0.2606    0.1503    0.0000 C   0  0
    4.9394    0.1503    0.0000 C   0  0
    3.6432   -2.0994    0.0000 C   0  0
   24.5681   -3.9951    0.0000 N   0  0
   25.8942   -3.2923    0.0000 N   0  3
   26.9545   -2.7304    0.0000 N   0  5
    8.7399    3.7207    0.0000 C   0  0
   10.2265    3.9945    0.0000 N   0  0
    8.0174    5.0180    0.0000 N   0  0
   10.4084    5.4994    0.0000 C   0  0  1  0  0  0
    9.0572    6.1153    0.0000 C   0  0  2  0  0  0
   10.6972    7.8929    0.0000 S   0  0
   11.4240    6.5733    0.0000 C   0  0  1  0  0  0
    9.2402    7.5918    0.0000 C   0  0
    8.2285    2.6351    0.0000 O   0  0
   12.9117    6.3752    0.0000 C   0  0
    8.1528    6.5421    0.0000 H   0  0
   11.3136    5.0744    0.0000 H   0  0
   13.4848    4.9881    0.0000 C   0  0
   14.9725    4.7900    0.0000 C   0  0
   15.5457    3.4029    0.0000 C   0  0
   17.0334    3.2048    0.0000 C   0  0
   17.6065    1.8177    0.0000 N   0  5
   17.7650    4.1560    0.0000 O   0  0
   19.0942    1.6196    0.0000 S   0  0
   19.6674    0.2326    0.0000 N   0  0
   18.6365    2.7289    0.0000 O   0  0
   19.8258    2.5708    0.0000 O   0  0
   21.1551    0.0345    0.0000 C   0  0
   21.7272   -1.3502    0.0000 C   0  0  2  0  0  0
   23.1848   -1.7047    0.0000 C   0  0  2  0  0  0
   23.2981   -3.2004    0.0000 C   0  0  1  0  0  0
   21.9197   -3.7468    0.0000 C   0  0  1  0  0  0
   20.9398   -2.6269    0.0000 O   0  0
   21.4586   -5.1750    0.0000 N   0  0
   21.4586   -7.5816    0.0000 N   0  0
   22.3337   -6.3601    0.0000 C   0  0
   20.0365   -7.1258    0.0000 C   0  0
   20.0365   -5.6308    0.0000 C   0  0
   18.7421   -4.8650    0.0000 N   0  0
   17.4294   -5.6308    0.0000 C   0  0
   17.4294   -7.1258    0.0000 N   0  0
   18.7421   -7.8915    0.0000 C   0  0
   18.7457   -9.0915    0.0000 N   0  0
   24.0995   -0.9280    0.0000 O   0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  2  5  1  0
  3  6  1  0
  4  7  1  0
  8  9  2  0
  9 10  2  0
 12 11  1  0
 13 11  1  0
 14 12  1  0
 15 13  1  0
 16 18  1  0
 17 14  1  0
 18 15  1  0
 19 11  2  0
 17 20  1  1
 15 21  1  6
 14 22  1  6
 15 14  1  0
 16 17  1  0
 20 23  1  0
 23 24  1  0
 24 25  1  0
 25 26  1  0
 26 27  1  0
 26 28  2  0
 27 29  1  0
 29 30  1  0
 29 31  2  0
 29 32  2  0
 30 33  1  0
 34 33  1  6
 34 35  1  0
 35 36  1  0
 36 37  1  0
 37 38  1  0
 38 34  1  0
 37 39  1  6
 39 43  1  0
 42 40  1  0
 40 41  2  0
 41 39  1  0
 42 43  2  0
 42 47  1  0
 43 44  1  0
 44 45  2  0
 45 46  1  0
 46 47  2  0
 47 48  1  0
 35 49  1  1
 36  8  1  6
M  CHG  4   1   1   9   1  10  -1  27  -1
M  END'''
        omolb = standardizer.standardize_molblock(molb)
        # print(omolb)
        self.assertEqual(omolb.find('1  2  2  3'), -1)

    def testValencePropagation2(self):
        ''' another example found in testing '''
        # first case:
        molb = '''
  SciTegic12231509382D

  3  2  0  0  0  0            999 V2000
    6.2708   -7.0500    0.0000 C   0  0
    7.0958   -7.0500    0.0000 N   0  0
    5.4458   -7.0500    0.0000 Se  0  5
  2  1  3  0
  3  1  1  0
M  CHG  1   3  -1
M  END'''
        omolb = standardizer.standardize_molblock(molb)
        # print(omolb)
        self.assertNotEqual(omolb.find('Se  0  0  0  0  0  2'), -1)
        self.assertEqual(Chem.MolBlockToInchi(omolb),
                         'InChI=1S/CHNSe/c2-1-3/h3H')

        molb = '''
  SciTegic11141416502D

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 S   0  0  0  0  0 15
M  END
'''
        omolb = standardizer.standardize_molblock(molb)
        self.assertEqual(Chem.MolBlockToInchi(omolb), 'InChI=1S/S')


#     def testUnchargingError(self):
#         ''' another example found in testing '''
#         # first case:
#         molb = '''
#   SciTegic11261410592D

#  14 13  0  0  0  0            999 V2000
#     0.0000    0.0000    0.0000 Br  0  0
#     2.1214   -0.3267    0.0000 O   0  5
#     2.5630    1.0323    0.0000 O   0  0
#     4.1866    0.4553    0.0000 N   0  3
#     3.7741   -0.8142    0.0000 N   0  3
#     2.7345    0.2253    0.0000 C   0  0
#     4.8540   -0.0296    0.0000 C   0  0
#     4.5991   -0.8142    0.0000 C   0  0
#     4.7386    1.0684    0.0000 C   0  0
#     3.6346    1.0684    0.0000 C   0  0
#     3.8603   -1.6347    0.0000 C   0  0
#     2.9671   -0.9858    0.0000 C   0  0
#     3.5191   -0.0885    0.0000 B   0  5
#     2.7419   -0.2461    0.0000 H   0  0
#   2  6  1  0
#   3  6  1  0
#   4  7  1  0
#   4  9  1  0
#   4 10  1  0
#   4 13  1  0
#   5  8  1  0
#   5 11  1  0
#   5 12  1  0
#   5 13  1  0
#   6 13  1  0
#   7  8  1  0
#  13 14  1  0
# M  CHG  4   2  -1   4   1   5   1  13  -1
# M  END'''
#         omolb = standardizer.standardize_molblock(molb)
#         print(omolb)
#         self.assertNotEqual(
#             omolb.find('M  CHG  4   2  -1   4   1   5   1  13  -1'), -1)
#         self.assertEqual(
#             Chem.MolBlockToInchi(omolb),
#             'InChI=1S/C7H19BN2O2.BrH/c1-9(2)5-6-10(3,4)8(9)7(11)12;/h7-8,11H,5-6H2,1-4H3;1H'
#         )

    def testGithub5(self):
        ''' standardizer.standardize_molblock crashes with some structures '''
        # first case:
        molb = '''


 37 36  0  0  0  0            999 V2000
   -0.7145   -2.8050    0.0000 C   0  0
    2.8579    0.9075    0.0000 C   0  0
   -0.7145   -1.9800    0.0000 C   0  0
    0.7145    0.4950    0.0000 C   0  0
   -0.7145    0.4950    0.0000 C   0  0
    0.0000    0.9075    0.0000 C   0  0
    0.7145   -0.3300    0.0000 C   0  0
   -0.7145   -0.3300    0.0000 C   0  0
    0.0000   -0.7425    0.0000 C   0  0
    1.4289    0.9075    0.0000 C   0  0
   -1.4289    0.9075    0.0000 C   0  0
    0.0000    1.7325    0.0000 I   0  0
    1.4289   -0.7425    0.0000 I   0  0
   -1.4289   -0.7425    0.0000 I   0  0
    2.1434    0.4950    0.0000 N   0  0
    0.0000   -1.5675    0.0000 N   0  0
   -1.4289   -1.5675    0.0000 O   0  0
    1.4289    1.7325    0.0000 O   0  0
   -2.1434    0.4950    0.0000 O   0  0
   -1.4289    1.7325    0.0000 O   0  0
   10.6951   -0.1941    0.0000 C   0  0
    9.2661   -0.1941    0.0000 C   0  0
    5.6938    0.2184    0.0000 C   0  0
    8.5517    0.2184    0.0000 C   0  0  1  0  0  0
    6.4083   -0.1941    0.0000 C   0  0  2  0  0  0
    7.8372   -0.1941    0.0000 C   0  0  2  0  0  0
    7.1227    0.2184    0.0000 C   0  0  2  0  0  0
    9.9806    0.2184    0.0000 N   0  0
    4.9793   -0.1941    0.0000 O   0  0
    8.5517    1.0434    0.0000 O   0  0
    6.4083   -1.0191    0.0000 O   0  0
    7.8372   -1.0191    0.0000 O   0  0
    7.1227    1.0434    0.0000 O   0  0
    8.5517   -0.6066    0.0000 H   0  0
    6.4083    0.6309    0.0000 H   0  0
    7.8372    0.6309    0.0000 H   0  0
    7.1227   -0.6066    0.0000 H   0  0
  3  1  1  0
  6  4  2  0
  6  5  1  0
  7  4  1  0
  8  5  2  0
  9  7  2  0
  9  8  1  0
 10  4  1  0
 11  5  1  0
 12  6  1  0
 13  7  1  0
 14  8  1  0
 15  2  1  0
 15 10  2  0
 16  3  2  0
 16  9  1  0
 17  3  1  0
 18 10  1  0
 19 11  2  0
 20 11  1  0
 24 22  1  0
 25 23  1  0
 26 24  1  0
 27 25  1  0
 27 26  1  0
 28 21  1  0
 28 22  1  0
 29 23  1  0
 24 30  1  6
 25 31  1  1
 26 32  1  6
 27 33  1  6
 24 34  1  6
 25 35  1  1
 26 36  1  6
 27 37  1  6
M  END'''
        omolb = standardizer.standardize_molblock(molb)
        self.assertEqual(
            Chem.MolBlockToInchi(omolb),
            'InChI=1S/C11H9I3N2O4.C7H17NO5/c1-3(17)16-9-7(13)4(10(18)15-2)6(12)5(8(9)14)11(19)20;1-8-2-4(10)6(12)7(13)5(11)3-9/h1-2H3,(H,15,18)(H,16,17)(H,19,20);4-13H,2-3H2,1H3/t;4-,5+,6+,7+/m.0/s1'
        )
        molb = """
  SciTegic07111214002D

 22 19  0  0  0  0            999 V2000
   12.9184   -8.1544    0.0000 C   0  0
   13.6343   -7.7410    0.0000 N   0  0
   12.2024   -7.7410    0.0000 O   0  0
   12.9184   -8.9811    0.0000 N   0  0
   14.3503   -8.1544    0.0000 C   0  0
   15.0663   -7.7410    0.0000 C   0  0
   15.7784   -8.1460    0.0000 C   0  0
   16.4924   -7.7327    0.0000 C   0  0
   17.2064   -8.1418    0.0000 C   0  0
   17.9204   -7.7285    0.0000 C   0  0
   18.6343   -8.1376    0.0000 C   0  0
   19.3483   -7.7243    0.0000 C   0  0
   20.0623   -8.1335    0.0000 C   0  0
   20.7763   -7.7202    0.0000 C   0  0
   21.4902   -8.1293    0.0000 C   0  0
   22.2042   -7.7160    0.0000 C   0  0
   22.9207   -8.1284    0.0000 N   0  0
   23.6362   -7.7143    0.0000 C   0  0
   24.3526   -8.1267    0.0000 O   0  0
   23.6352   -6.8875    0.0000 N   0  0
   14.3464   -5.9039    0.0000 Cl  0  0
   24.1332   -6.2672    0.0000 Cl  0  0
  2  5  1  0
 11 12  1  0
  1  2  1  0
 12 13  1  0
  5  6  1  0
 13 14  1  0
 14 15  1  0
  6  7  1  0
 15 16  1  0
  1  3  1  0
 16 17  1  0
  7  8  1  0
 17 18  1  0
 18 19  1  0
  8  9  1  0
 18 20  2  0
  1  4  2  0
  9 10  1  0
 10 11  1  0
M  END"""
        omolb = standardizer.standardize_molblock(molb)
        self.assertEqual(
            Chem.MolBlockToInchi(omolb),
            'InChI=1S/C14H30N4O2.2ClH/c15-13(19)17-11-9-7-5-3-1-2-4-6-8-10-12-18-14(16)20;;/h1-12H2,(H3,15,17,19)(H3,16,18,20);2*1H'
        )

    def testGithub6(self):
        ''' standardiser changes the structure even with excluded structures '''
        # first case:
        molb = '''ranitidine bismuth citrate
  SciTegic03151315062D

 35 33  0  0  0  0            999 V2000
    5.1714   -6.6800    0.0000 C   0  0
    4.4533   -7.0848    0.0000 C   0  0
    5.8848   -7.0848    0.0000 C   0  0
    4.7570   -5.9667    0.0000 C   0  0
    6.5932   -6.6800    0.0000 C   0  0
    4.4678   -7.8945    0.0000 C   0  0
    5.1570   -5.2486    0.0000 O   0  0
    5.2003   -8.2849    0.0000 O   0  0
    6.5932   -5.8558    0.0000 O   0  0
    5.5763   -5.9667    0.0000 O   0  0
    3.9328   -5.9715    0.0000 O   0  5
    3.7738   -8.3283    0.0000 O   0  5
    7.3113   -7.0848    0.0000 O   0  5
    7.1571   -8.3283    0.0000 Bi  0  1
    2.0387   -7.2005    0.0000 N   0  3
    1.3591   -6.7812    0.0000 C   0  0
    1.3591   -6.0149    0.0000 C   0  0
   -4.2943   -6.0389    0.0000 C   0  0
   -3.6292   -5.5956    0.0000 O   0  0
   -2.9303   -6.0389    0.0000 C   0  0
   -4.0340   -6.8535    0.0000 C   0  0
   -3.1905   -6.8535    0.0000 C   0  0
    2.7472   -6.8197    0.0000 O   0  5
    2.0580   -7.9813    0.0000 O   0  0
   -5.0172   -5.5956    0.0000 C   0  0
    2.0725   -5.5956    0.0000 N   0  0
    0.6121   -5.5956    0.0000 N   0  0
   -5.7257   -6.0389    0.0000 N   0  0
   -1.5037   -6.0149    0.0000 S   0  0
   -2.2218   -5.5956    0.0000 C   0  0
   -0.0819   -6.0149    0.0000 C   0  0
   -0.8097   -5.5956    0.0000 C   0  0
   -5.7257   -6.8535    0.0000 C   0  0
   -6.4438   -5.5956    0.0000 C   0  0
    2.7665   -6.0149    0.0000 C   0  0
  6  2  1  0
  7  4  2  0
  8  6  2  0
  9  5  2  0
 10  1  1  0
 11  4  1  0
 12  6  1  0
 13  5  1  0
  2  1  1  0
  3  1  1  0
  4  1  1  0
  5  3  1  0
 16 15  1  0
 17 16  2  3
 18 19  1  0
 19 20  1  0
 20 30  1  0
 21 22  1  0
 22 20  2  0
 23 15  1  0
 24 15  2  0
 25 18  1  0
 26 17  1  0
 27 17  1  0
 28 25  1  0
 29 32  1  0
 30 29  1  0
 31 27  1  0
 32 31  1  0
 33 28  1  0
 34 28  1  0
 35 26  1  0
 21 18  2  0
M  CHG  6  11  -1  12  -1  13  -1  14   3  15   1  23  -1
M  END'''
        m = Chem.MolFromMolBlock(molb, sanitize=False, removeHs=False)
        exclude = standardizer.exclude_flag(m, includeRDKitSanitization=False)
        self.assertTrue(exclude)
        omolb = standardizer.standardize_molblock(molb)
        inchi = Chem.MolBlockToInchi(molb)
        self.assertEqual(Chem.MolBlockToInchi(omolb), inchi)

        # check salt stripping for this
        m = Chem.MolFromMolBlock(molb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 3)
        smb, exclude = standardizer.get_parent_molblock(molb)
        sm = Chem.MolFromMolBlock(smb)
        print(Chem.MolToSmiles(sm))
        self.assertEqual(len(Chem.GetMolFrags(sm)), 1)

    def testDoubleBondStereoAndTransforms(self):
        molb = '''
  Mrv1810 10111911132D          

  7  6  0  0  0  0            999 V2000
    0.2030   -1.7143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9175   -1.3018    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5115   -1.3018    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.2030   -2.5393    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2259   -1.7143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9404   -1.3018    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6549   -1.7143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  2  0  0  0  0
  3  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  7  1  0  0  0  0
M  END
'''
        smb = standardizer.standardize_molblock(molb)
        self.assertEqual(smb.find('5  6  2  3'), -1)

    def testDuplicatesPostStripping(self):
        molb = '''duplicates, no chirality
  Mrv1824 10241910462D          

  9  6  0  0  0  0            999 V2000
   -3.4152    1.6964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7007    2.1089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9862    1.6964    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -2.7007    2.9339    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.4239    1.5964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1384    2.0089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8529    1.5964    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    1.1384    2.8339    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9152    2.6563    0.0000 Ca  0  2  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  6  8  1  0  0  0  0
M  CHG  3   3  -1   7  -1   9   2
M  END
'''
        smb, exclude = standardizer.get_parent_molblock(molb)
        m = Chem.MolFromMolBlock(smb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 1)

        molb = '''duplicates, chirality
  Mrv1824 10241910462D          

  9  6  0  0  0  0            999 V2000
   -3.4152    1.6964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7007    2.1089    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -1.9862    1.6964    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -2.7007    2.9339    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.4239    1.5964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1384    2.0089    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    1.8529    1.5964    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    1.1384    2.8339    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9152    2.6563    0.0000 Ca  0  2  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  2  4  1  1  0  0  0
  6  8  1  1  0  0  0
M  CHG  3   3  -1   7  -1   9   2
M  END
'''
        smb, exclude = standardizer.get_parent_molblock(molb)
        m = Chem.MolFromMolBlock(smb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 1)
        molb = '''duplicates except for chirality
  Mrv1824 10241910462D          

  9  6  0  0  0  0            999 V2000
   -3.4152    1.6964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7007    2.1089    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -1.9862    1.6964    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -2.7007    2.9339    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.4239    1.5964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1384    2.0089    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    1.8529    1.5964    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    1.1384    2.8339    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9152    2.6563    0.0000 Ca  0  2  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  2  4  1  1  0  0  0
  6  8  1  6  0  0  0
M  CHG  3   3  -1   7  -1   9   2
M  END
'''
        smb, exclude = standardizer.get_parent_molblock(molb)
        m = Chem.MolFromMolBlock(smb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 2)

    def testCovalentAlkalineNBonds(self):
        mb = '''
  Mrv1810 10301916122D          

  6  6  0  0  0  0            999 V2000
   -0.3835    0.3550    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0979    0.7675    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0979    1.5925    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3310    1.5925    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.3310    0.7675    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3835   -0.4700    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  4  5  2  0  0  0  0
  1  5  1  0  0  0  0
  3  4  1  0  0  0  0
  1  6  1  0  0  0  0
M  END
'''
        smb = standardizer.standardize_molblock(mb)
        m = Chem.MolFromMolBlock(smb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 2)
        self.assertEqual(m.GetAtomWithIdx(0).GetAtomicNum(), 11)
        self.assertEqual(m.GetAtomWithIdx(0).GetFormalCharge(), 1)
        self.assertEqual(m.GetAtomWithIdx(1).GetAtomicNum(), 7)
        self.assertEqual(m.GetAtomWithIdx(1).GetFormalCharge(), -1)

    def testGetParentOrder(self):
        ' make sure we remove isotopes before stripping salts '
        mb = '''
  Mrv1810 11071914272D          

 10  9  0  0  0  0            999 V2000
   -5.3684    3.1483    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0829    2.7358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0829    1.9108    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3684    1.4983    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6540    1.9108    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6540    2.7358    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3684    3.9733    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5188    3.8797    0.0000 D   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8043    4.2922    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0899    3.8797    0.0000 D   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  1  6  2  0  0  0  0
  1  7  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  1  0  0  0  0
M  END
'''
        m = Chem.MolFromMolBlock(mb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 2)
        smb, exclude = standardizer.get_parent_molblock(mb)
        m = Chem.MolFromMolBlock(smb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 1)

    def testGetChargedSalts(self):
        ' duplicates only after uncharging '
        mb = '''1249071
  SciTegic10031215492D

 94 94  0  0  0  0            999 V2000
   -0.9320  -13.0513    0.0000 C   0  0
   -1.6459  -12.6417    0.0000 C   0  0
   -1.6481  -11.8185    0.0000 C   0  0
   -0.9364  -11.4050    0.0000 C   0  0
   -0.2223  -11.8147    0.0000 C   0  0
   -0.2201  -12.6378    0.0000 C   0  0
    0.4938  -13.0475    0.0000 C   0  0
    1.2056  -12.6340    0.0000 N   0  0
    0.4961  -13.8706    0.0000 O   0  0
    1.9197  -13.0436    0.0000 C   0  0
    2.6314  -12.6301    0.0000 C   0  0
    3.3454  -13.0398    0.0000 N   0  0
    2.6291  -11.8070    0.0000 O   0  0
    4.0572  -12.6263    0.0000 C   0  0
    4.0550  -11.8031    0.0000 C   0  0
    4.7667  -11.3896    0.0000 C   0  0
    5.4807  -11.7993    0.0000 C   0  0
    5.4829  -12.6225    0.0000 C   0  0
    4.7712  -13.0360    0.0000 C   0  0
    4.7734  -13.8591    0.0000 I   0  0
    3.3409  -11.3935    0.0000 I   0  0
    6.1925  -11.3858    0.0000 I   0  0
    4.7645  -10.5665    0.0000 C   0  0
    5.4763  -10.1530    0.0000 N   0  0
    4.0505  -10.1568    0.0000 O   0  0
    6.1903  -10.5626    0.0000 C   0  0
    6.9021  -10.1491    0.0000 C   0  0
    7.6160  -10.5588    0.0000 O   0  5
   -0.9386  -10.5818    0.0000 C   0  0
   -0.2268  -10.1683    0.0000 O   0  0
   -1.6525  -10.1722    0.0000 N   0  0
   -2.3643  -10.5857    0.0000 C   0  0
    0.4894  -11.4011    0.0000 I   0  0
   -2.3621  -11.4088    0.0000 I   0  0
   -2.3577  -13.0552    0.0000 N   0  0
   -2.3555  -13.8783    0.0000 C   0  0
   -3.0717  -12.6455    0.0000 C   0  0
   -3.7835  -13.0590    0.0000 C   0  0
   -3.0739  -11.8223    0.0000 O   0  0
   -0.9298  -13.8745    0.0000 I   0  0
   -0.5948   -5.8212    0.0000 C   0  0
   -1.3089   -5.4115    0.0000 C   0  0
   -1.3111   -4.5884    0.0000 C   0  0
   -0.5993   -4.1749    0.0000 C   0  0
    0.1147   -4.5845    0.0000 C   0  0
    0.1170   -5.4077    0.0000 C   0  0
    0.8310   -5.8173    0.0000 C   0  0
    1.5428   -5.4038    0.0000 N   0  0
    0.8332   -6.6405    0.0000 O   0  0
    2.2567   -5.8135    0.0000 C   0  0
    2.9685   -5.4000    0.0000 C   0  0
    3.6825   -5.8097    0.0000 N   0  0
    2.9663   -4.5768    0.0000 O   0  0
    4.3942   -5.3962    0.0000 C   0  0
    4.3920   -4.5730    0.0000 C   0  0
    5.1038   -4.1595    0.0000 C   0  0
    5.8178   -4.5692    0.0000 C   0  0
    5.8200   -5.3924    0.0000 C   0  0
    5.1082   -5.8058    0.0000 C   0  0
    5.1104   -6.6290    0.0000 I   0  0
    3.6781   -4.1633    0.0000 I   0  0
    6.5296   -4.1557    0.0000 I   0  0
    5.1016   -3.3363    0.0000 C   0  0
    5.8134   -2.9228    0.0000 N   0  0
    4.3877   -2.9267    0.0000 O   0  0
    6.5274   -3.3325    0.0000 C   0  0
    7.2391   -2.9190    0.0000 C   0  0
    7.9531   -3.3286    0.0000 O   0  0
   -0.6015   -3.3517    0.0000 C   0  0
    0.1103   -2.9382    0.0000 O   0  0
   -1.3154   -2.9420    0.0000 N   0  0
   -2.0272   -3.3555    0.0000 C   0  0
    0.8265   -4.1710    0.0000 I   0  0
   -2.0250   -4.1787    0.0000 I   0  0
   -2.0206   -5.8250    0.0000 N   0  0
   -2.0183   -6.6482    0.0000 C   0  0
   -2.7346   -5.4154    0.0000 C   0  0
   -3.4464   -5.8289    0.0000 C   0  0
   -2.7368   -4.5922    0.0000 O   0  0
   -0.5926   -6.6443    0.0000 I   0  0
    7.4443  -12.1174    0.0000 Na  0  3
    9.2524   -6.7030    0.0000 N   0  0
    9.9664   -7.1127    0.0000 C   0  0
   13.5297   -6.6915    0.0000 C   0  0
   14.2437   -7.1012    0.0000 O   0  0
    8.5406   -7.1165    0.0000 C   0  0
   10.6759   -5.8760    0.0000 O   0  0
   12.1017   -5.8722    0.0000 O   0  0
   11.3944   -7.9320    0.0000 O   0  0
   12.8201   -7.9282    0.0000 O   0  0
   10.6792   -6.7000    0.0000 C   0  0  2  0  0  0
   11.3917   -7.1083    0.0000 C   0  0  2  0  0  0
   12.1042   -6.6958    0.0000 C   0  0  1  0  0  0
   12.8167   -7.1042    0.0000 C   0  0  1  0  0  0
  1  2  1  0
  1  6  2  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  1  0
  7  8  1  0
  7  9  2  0
  8 10  1  0
 10 11  1  0
 11 12  1  0
 11 13  2  0
 12 14  1  0
 14 15  1  0
 14 19  2  0
 15 16  2  0
 16 17  1  0
 17 18  2  0
 18 19  1  0
 19 20  1  0
 15 21  1  0
 17 22  1  0
 16 23  1  0
 23 24  1  0
 23 25  2  0
 24 26  1  0
 26 27  1  0
 27 28  1  0
  4 29  1  0
 29 30  2  0
 29 31  1  0
 31 32  1  0
  5 33  1  0
  3 34  1  0
  2 35  1  0
 35 36  1  0
 35 37  1  0
 37 38  1  0
 37 39  2  0
  1 40  1  0
 41 42  1  0
 41 46  2  0
 42 43  2  0
 43 44  1  0
 44 45  2  0
 45 46  1  0
 46 47  1  0
 47 48  1  0
 47 49  2  0
 48 50  1  0
 50 51  1  0
 51 52  1  0
 51 53  2  0
 52 54  1  0
 54 55  1  0
 54 59  2  0
 55 56  2  0
 56 57  1  0
 57 58  2  0
 58 59  1  0
 59 60  1  0
 55 61  1  0
 57 62  1  0
 56 63  1  0
 63 64  1  0
 63 65  2  0
 64 66  1  0
 66 67  1  0
 67 68  1  0
 44 69  1  0
 69 70  2  0
 69 71  1  0
 71 72  1  0
 45 73  1  0
 43 74  1  0
 42 75  1  0
 75 76  1  0
 75 77  1  0
 77 78  1  0
 77 79  2  0
 41 80  1  0
 82 83  1  0
 84 85  1  0
 82 86  1  0
 94 84  1  0
 91 92  1  0
 91 87  1  6
 92 93  1  0
 92 89  1  6
 83 91  1  0
 93 88  1  6
 93 94  1  0
 94 90  1  1
M  CHG  2  28  -1  81   1
M  END'''
        m = Chem.MolFromMolBlock(mb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 4)
        smb, exclude = standardizer.get_parent_molblock(mb)
        m = Chem.MolFromMolBlock(smb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 1)

    def testExcludeFlagAndSalts(self):
        ' exclude flag and salts '
        mb = '''674781
  SciTegic01161211252D

  6  4  0  0  0  0            999 V2000
    5.6417   -7.4167    0.0000 S   0  0
    5.6417   -8.1292    0.0000 O   0  5
    4.9375   -7.4167    0.0000 O   0  5
    6.3542   -7.4167    0.0000 O   0  0
    5.6417   -6.7125    0.0000 O   0  0
    4.2125   -8.3750    0.0000 Fe  0  2
  2  1  1  0
  3  1  1  0
  4  1  2  0
  5  1  2  0
M  CHG  3   2  -1   3  -1   6   2
M  END'''
        m = Chem.MolFromMolBlock(mb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 2)
        smb, exclude = standardizer.get_parent_molblock(mb)
        sm = Chem.MolFromMolBlock(smb)
        self.assertEqual(len(Chem.GetMolFrags(sm)), 2)
        self.assertEqual(sm.GetNumAtoms(), m.GetNumAtoms())
        self.assertTrue(exclude)
        smb, exclude = standardizer.get_parent_molblock(mb,
                                                        check_exclusion=False)
        sm = Chem.MolFromMolBlock(smb)
        self.assertEqual(len(Chem.GetMolFrags(sm)), 1)
        self.assertEqual(sm.GetNumAtoms(), 1)
        self.assertFalse(exclude)

        mb = '''449259
  SciTegic12011715492D

  5  0  0  0  0  0            999 V2000
    0.0295    1.3413    0.0000 O   0  0
   -0.6707    0.0147    0.0000 O   0  0
    1.3561    0.6412    0.0000 Cl  0  5
    0.6559   -0.6854    0.0000 Cl  0  5
   -1.3708   -1.3118    0.0000 Cu  0  2
M  CHG  3   3  -1   4  -1   5   2
M  END'''
        m = Chem.MolFromMolBlock(mb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 5)
        smb, exclude = standardizer.get_parent_molblock(mb)
        sm = Chem.MolFromMolBlock(smb, sanitize=False)
        self.assertEqual(len(Chem.GetMolFrags(sm)), 3)
        self.assertTrue(exclude)

        mb = '''453882
  SciTegic12011715492D

 10  0  0  0  0  0            999 V2000
    2.2357   -0.2536    0.0000 Fe  0  1
    1.1448    0.7759    0.0000 Cl  0  5
    0.0538    1.8054    0.0000 Cl  0  5
    1.2061   -1.3445    0.0000 Cl  0  5
    0.1152   -0.3150    0.0000 O   0  0
   -0.9757    0.7145    0.0000 O   0  0
    0.1766   -2.4354    0.0000 O   0  0
   -0.9143   -1.4059    0.0000 O   0  0
   -2.0052   -0.3764    0.0000 O   0  0
   -1.0371    2.8350    0.0000 O   0  0
M  CHG  4   1   3   2  -1   3  -1   4  -1
M  END'''
        m = Chem.MolFromMolBlock(mb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 10)
        smb, exclude = standardizer.get_parent_molblock(mb)
        sm = Chem.MolFromMolBlock(smb, sanitize=False)
        self.assertEqual(len(Chem.GetMolFrags(sm)), 4)
        self.assertTrue(exclude)

    def testChiralHsAndSalts(self):
        ' check that salts with chiral Hs are removed '
        mb = '''674457
  Mrv0541 01111613452D          

 80 82  0  0  0  0            999 V2000
    5.7088  -16.3933    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    1.0826  -16.8108    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -1.5534  -15.5293    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -2.9806  -15.5389    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3010  -17.1021    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -3.6942  -13.4855    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -2.2670  -15.9467    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    0.6748  -17.5244    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    4.4758  -17.1021    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -0.8398  -14.3059    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -1.5534  -14.7137    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    7.1409  -16.3787    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    5.3010  -15.6797    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.8496  -15.9564    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -2.9806  -12.2622    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    8.5681  -16.3641    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -2.9806  -14.7137    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9806  -13.0826    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    6.4224  -15.9758    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2670  -11.8448    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0777  -15.7137    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0582  -16.3933    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -3.1699  -19.1264    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -0.1553  -17.5244    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -0.8398  -13.4855    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -1.9369  -18.4128    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    4.4758  -15.6797    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -2.3447  -19.1264    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6748  -16.0972    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3447  -17.6992    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6942  -14.3059    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -3.5874  -18.4128    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    9.2768  -15.9418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8398  -15.9418    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1699  -17.6992    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -1.5534  -12.2525    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -1.5534  -13.0777    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5680  -16.8108    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5631  -18.2380    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    7.1456  -17.2040    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -3.6942  -15.9516    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1553  -16.0972    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -2.2719  -11.0243    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.2719  -15.1166    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1068  -12.7622    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7137  -17.8205    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0826  -18.2380    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0631  -17.8205    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.5777  -17.1895    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.8399  -15.1312    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4320  -12.7622    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2330  -16.3933    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2670  -13.4952    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8884  -19.5341    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4126  -18.4128    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.9952  -16.3496    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5146  -13.4855    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2622  -16.7671    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4369  -17.6215    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2719  -14.2962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6942  -11.8495    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0631  -14.9613    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1748  -19.9516    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0146  -13.4855    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4078  -14.7234    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2379  -14.9613    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.8738  -18.4321    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5825  -16.9807    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8398  -11.8399    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.8641  -17.6069    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1505  -18.9516    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3884  -18.2380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5631  -15.3787    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6020  -19.1118    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1214  -14.3107    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7962  -17.2186    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0146  -14.3059    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2233  -18.8205    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    7.8545  -16.7914    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.4224  -16.8059    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  3 11  1  0  0  0  0
  4  7  1  0  0  0  0
  5  1  1  0  0  0  0
  6 31  1  0  0  0  0
  7  3  1  0  0  0  0
  8  2  1  0  0  0  0
  9  5  1  0  0  0  0
 10 21  1  0  0  0  0
 11 10  1  0  0  0  0
 12 19  1  0  0  0  0
 13  1  1  0  0  0  0
 14 12  1  0  0  0  0
 15 20  1  0  0  0  0
 16 14  1  0  0  0  0
 17  4  1  0  0  0  0
 18 15  1  0  0  0  0
 19  1  1  0  0  0  0
 20 36  1  0  0  0  0
 21  2  1  0  0  0  0
 22 27  1  0  0  0  0
 23 28  1  0  0  0  0
 24  8  1  0  0  0  0
 25 10  1  0  0  0  0
 26 34  1  0  0  0  0
 27 13  1  0  0  0  0
 28 26  1  0  0  0  0
 29  2  1  0  0  0  0
 30 26  1  0  0  0  0
 31 17  1  0  0  0  0
 32 35  1  0  0  0  0
 33 16  1  0  0  0  0
  3 34  1  6  0  0  0
 35 30  1  0  0  0  0
 36 37  1  0  0  0  0
 37 25  1  0  0  0  0
 38 42  1  0  0  0  0
 24 39  1  1  0  0  0
 40 12  1  0  0  0  0
 41  4  2  0  0  0  0
 42 29  1  0  0  0  0
 43 20  2  0  0  0  0
 44 33  2  0  0  0  0
  6 45  1  1  0  0  0
  5 46  1  6  0  0  0
  8 47  1  6  0  0  0
  9 48  1  1  0  0  0
 16 49  1  6  0  0  0
 14 50  1  6  0  0  0
 25 51  1  1  0  0  0
 22 52  1  1  0  0  0
 18 53  1  1  0  0  0
 23 54  1  1  0  0  0
 32 55  1  1  0  0  0
 56 33  1  0  0  0  0
 57  6  1  0  0  0  0
  7 58  1  1  0  0  0
 40 59  1  1  0  0  0
 11 60  1  6  0  0  0
 15 61  1  1  0  0  0
 27 62  1  1  0  0  0
 63 23  1  0  0  0  0
 64 25  1  0  0  0  0
 31 65  1  6  0  0  0
 66 62  1  0  0  0  0
 67 70  1  0  0  0  0
 35 68  1  6  0  0  0
 36 69  1  6  0  0  0
 70 40  1  0  0  0  0
 71 39  1  0  0  0  0
 72 39  1  0  0  0  0
 42 73  1  1  0  0  0
 74 54  1  0  0  0  0
 75 65  1  0  0  0  0
  2 76  1  6  0  0  0
 10 77  1  1  0  0  0
 26 78  1  6  0  0  0
 12 79  1  6  0  0  0
  1 80  1  6  0  0  0
 22  9  1  0  0  0  0
 38 24  1  0  0  0  0
  6 18  1  0  0  0  0
 32 23  1  0  0  0  0
M  END'''
        m = Chem.MolFromMolBlock(mb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 2)
        smb, exclude = standardizer.get_parent_molblock(mb)
        sm = Chem.MolFromMolBlock(smb)
        self.assertEqual(len(Chem.GetMolFrags(sm)), 1)

    def testSaltStrippingAromaticity(self):
        ' check that get_parent_molblock does not set aromaticity '
        mb = '''886596
  SciTegic03291108532D

 25 26  0  0  0  0            999 V2000
   -0.3248    1.8375    0.0000 O   0  0
   -1.7537   -0.6375    0.0000 O   0  0
    1.8187   -0.2250    0.0000 O   0  0
    1.1042   -1.4625    0.0000 O   0  0
    6.2014    0.7366    0.0000 N   0  0
    4.7724    0.7366    0.0000 N   0  0
    0.3897    0.6000    0.0000 N   0  0
   -1.0392    0.6000    0.0000 N   0  0
    6.2014   -1.7384    0.0000 N   0  0
    6.9158   -0.5009    0.0000 C   0  0
    6.9158    0.3241    0.0000 C   0  0
    6.2014   -0.9134    0.0000 C   0  0
    5.4869    0.3241    0.0000 C   0  0
    5.4869   -0.5009    0.0000 C   0  0
    0.3897   -0.2250    0.0000 C   0  0
    7.6303   -0.9134    0.0000 C   0  0
    7.6303    0.7366    0.0000 C   0  0
   -0.3248    1.0125    0.0000 C   0  0
   -1.0392   -0.2250    0.0000 C   0  0
   -0.3248   -0.6375    0.0000 C   0  0
    1.1042   -0.6375    0.0000 C   0  0
    8.3448   -0.5009    0.0000 C   0  0
    8.3448    0.3241    0.0000 C   0  0
    4.0579    0.3241    0.0000 C   0  0
    4.7724    1.5616    0.0000 C   0  0
  1 18  2  0
  2 19  2  0
  3 21  2  0
  4 21  1  0
  5 11  2  0
  5 13  1  0
  6 13  1  0
  6 24  1  0
  6 25  1  0
  7 15  1  0
  7 18  1  0
  8 18  1  0
  8 19  1  0
  9 12  1  0
 10 11  1  0
 10 12  2  0
 10 16  1  0
 11 17  1  0
 12 14  1  0
 13 14  2  0
 15 20  2  0
 15 21  1  0
 16 22  2  0
 17 23  2  0
 19 20  1  0
 22 23  1  0
M  END'''
        m = Chem.MolFromMolBlock(mb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 2)
        smb, exclude = standardizer.get_parent_molblock(mb)
        sm = Chem.MolFromMolBlock(smb, sanitize=False)
        self.assertEqual(len(Chem.GetMolFrags(sm)), 1)
        for bond in sm.GetBonds():
            self.assertFalse(bond.GetIsAromatic())

    def testSaltStrippingNumBonds(self):
        ' make sure the salt stripper actually matches the entire fragment, not just the atoms '
        mb = '''205714
          11280715242D 1   1.00000     0.00000     0

  9  8  0     0  0            999 V2000
    5.7250   -7.3292    0.0000 N   0  3  3  0  0  0           0  0  0
    6.1417   -8.0417    0.0000 C   0  0  0  0  0  0           0  0  0
    5.3125   -8.0417    0.0000 C   0  0  0  0  0  0           0  0  0
    7.8042   -7.0875    0.0000 Cl  0  5  0  0  0  0           0  0  0
    6.5042   -6.9792    0.0000 C   0  0  0  0  0  0           0  0  0
    5.0042   -6.9167    0.0000 C   0  0  0  0  0  0           0  0  0
    7.2125   -5.7417    0.0000 O   0  0  0  0  0  0           0  0  0
    6.5042   -6.1542    0.0000 C   0  0  0  0  0  0           0  0  0
    5.0042   -6.0917    0.0000 C   0  0  0  0  0  0           0  0  0
  2  1  1  0     0  0
  3  1  1  0     0  0
  5  1  1  0     0  0
  6  1  1  0     0  0
  7  8  1  0     0  0
  8  5  1  0     0  0
  9  6  1  0     0  0
  3  2  1  0     0  0
M  CHG  2   1   1   4  -1
M  END'''
        m = Chem.MolFromMolBlock(mb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 2)
        smb, exclude = standardizer.get_parent_molblock(mb)
        sm = Chem.MolFromMolBlock(smb, sanitize=False)
        self.assertEqual(len(Chem.GetMolFrags(sm)), 1)

    def testSolventSaltRemovalProblem(self):
        ' make sure the solvents are removed from mols where all frags are either solvents or salts '
        mb = '''1673174
 
 
 12 11  0     0  0  0  0  0  0999 V2000
    3.3001   -0.8250    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    3.3001   -1.6501    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.1251   -0.8250    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.3001    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8934   -2.5410    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4751   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0667   -1.5386    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0667   -0.1114    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8250   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2417   -1.5386    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2417   -0.1114    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  3  2  0  0  0  0
  1  4  1  0  0  0  0
  1  6  1  0  0  0  0
  6  7  2  0  0  0  0
  6  8  1  0  0  0  0
  7 10  1  0  0  0  0
  8 11  2  0  0  0  0
  9 10  2  0  0  0  0
  9 11  1  0  0  0  0
  9 12  1  0  0  0  0
M  END '''
        m = Chem.MolFromMolBlock(mb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 2)
        smb, exclude = standardizer.get_parent_molblock(mb)
        sm = Chem.MolFromMolBlock(smb, sanitize=False)
        self.assertEqual(len(Chem.GetMolFrags(sm)), 1)

    def testSaltRemovalAndCharges(self):
        ' make sure the whole remaining molecule is neutralized, not the fragments '
        mb = '''1674265
  SciTegic12231509382D

 17 16  0  0  0  0            999 V2000
    4.4613   -2.5540    0.0000 P   0  0
    4.4613   -1.7375    0.0000 O   0  0
    5.2844   -2.5540    0.0000 O   0  0
    4.4613   -3.3836    0.0000 O   0  5
    3.6383   -2.5540    0.0000 O   0  0
    1.4305   -1.0778    0.0000 N   0  0
    2.6977   -0.6598    0.0000 N   0  3
    2.2078   -1.3260    0.0000 C   0  0
    1.4305   -0.2482    0.0000 C   0  0
    2.2143    0.0000    0.0000 C   0  0
    5.7024   -3.2725    0.0000 C   0  0
    5.1734   -1.3260    0.0000 C   0  0
    0.7577   -1.5676    0.0000 C   0  0
    3.5273   -0.6598    0.0000 C   0  0
    6.5320   -3.2725    0.0000 C   0  0
    5.1734   -0.5095    0.0000 C   0  0
    0.0000   -1.2346    0.0000 C   0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  2  0
  2 12  1  0
  3 11  1  0
  6  8  1  0
  6  9  1  0
  6 13  1  0
  7  8  2  0
  7 10  1  0
  7 14  1  0
  9 10  2  0
 11 15  1  0
 12 16  1  0
 13 17  1  0
M  CHG  2   4  -1   7   1
M  END'''
        m = Chem.MolFromMolBlock(mb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 2)
        smb, exclude = standardizer.get_parent_molblock(mb)
        sm = Chem.MolFromMolBlock(smb, sanitize=False)
        self.assertEqual(len(Chem.GetMolFrags(sm)), 2)
        self.assertEqual(Chem.GetFormalCharge(sm), 0)

    def testLosingCrossedDoubleBond(self):
        ' cases where we lost a crossed double bond on standardization '
        mb = '''1239046
  SciTegic02161214032D

 24 28  0  0  0  0            999 V2000
    0.9032    3.5258    0.0000 Br  0  0
    2.7902    1.4518    0.0000 O   0  0
    3.6834    0.4598    0.0000 O   0  0
   -2.0248    2.0333    0.0000 O   0  0
   -2.0248    3.3682    0.0000 O   0  0
   -0.5488    0.8208    0.0000 O   0  0
    0.4907   -0.2188    0.0000 N   0  0
    1.5707    0.5658    0.0000 C   0  0
    0.9032    1.0508    0.0000 C   0  0
    2.3777    0.7374    0.0000 C   0  0
    1.3157   -0.2188    0.0000 C   0  0
    0.2358    0.5658    0.0000 C   0  0
    2.9297    0.1243    0.0000 C   0  0
    0.9032    1.8758    0.0000 C   0  0
    0.1888    2.2883    0.0000 C   0  0
    1.8678   -0.8319    0.0000 C   0  0
   -1.2402    2.2883    0.0000 C   0  0
    2.6747   -0.6604    0.0000 C   0  0
    0.1888    3.1133    0.0000 C   0  0
   -1.2402    3.1133    0.0000 C   0  0
   -0.5257    1.8758    0.0000 C   0  0
   -0.5257    3.5258    0.0000 C   0  0
    3.5971    1.2803    0.0000 C   0  0
   -2.5097    2.7008    0.0000 C   0  0
  1 19  1  0
  2 10  1  0
  2 23  1  0
  3 13  1  0
  3 23  1  0
  4 17  1  0
  4 24  1  0
  5 20  1  0
  5 24  1  0
  6 12  2  0
  7 11  1  0
  7 12  1  0
  8  9  1  0
  8 10  2  0
  8 11  1  0
  9 12  1  0
  9 14  2  3
 10 13  1  0
 11 16  2  0
 13 18  2  0
 14 15  1  0
 15 19  1  0
 15 21  2  0
 16 18  1  0
 17 20  2  0
 17 21  1  0
 19 22  2  0
 20 22  1  0
M  END'''
        m = Chem.MolFromMolBlock(mb)
        self.assertEqual(
            m.GetBondBetweenAtoms(8, 13).GetStereo(),
            Chem.BondStereo.STEREOANY)
        smb = standardizer.standardize_molblock(mb)
        sm = Chem.MolFromMolBlock(smb, sanitize=False)
        self.assertEqual(
            sm.GetBondBetweenAtoms(8, 13).GetStereo(),
            Chem.BondStereo.STEREOANY)


    def testChiralHsInterferingWithSaltStripping(self):
        ' cases where chiral Hs led to duplicates in salt stripper results '
        mb = '''codeine sulfate
  Mrv17a0 10081811232D
          
 57 61  0  0  0  0            999 V2000
   17.3075   -9.5099    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   16.4799   -9.5199    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   17.7189   -8.8027    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   16.0685  -10.2272    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   17.7038  -10.2272    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   16.8561  -10.6485    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   17.3126   -8.0954    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   16.0685   -8.8027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   18.5364   -8.8077    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   16.4799   -8.0903    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   17.7189   -7.3781    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   16.8862   -8.9281    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   18.9478   -9.5250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   18.5313  -10.2372    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   15.2459  -10.2272    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   16.8962   -7.3731    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   15.2459   -8.8027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   14.8346   -9.5199    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   18.9427  -10.9596    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   14.8246  -10.9446    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   18.4461   -6.9617    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   14.0019  -10.9446    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   18.3251   -7.7627    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   21.7639  -10.4057    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   21.7639   -9.5529    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   21.7639  -11.2685    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   20.9011  -10.4057    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   22.6115  -10.4057    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   16.7096   -7.0535    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   17.9719  -11.4008    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   17.4616  -14.5725    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   16.6340  -14.5825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   17.8730  -13.8652    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   16.2226  -15.2898    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   17.8579  -15.2898    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   17.0102  -15.7111    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   17.4667  -13.1580    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   16.2226  -13.8652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   18.6905  -13.8703    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   16.6340  -13.1529    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   17.8730  -12.4407    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   17.0403  -13.9907    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   19.1019  -14.5875    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   18.6855  -15.2997    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   15.4000  -15.2898    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   17.0503  -12.4357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   15.4000  -13.8652    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   14.9887  -14.5825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   19.0968  -16.0222    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   14.9787  -16.0071    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   18.6002  -12.0243    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   14.1560  -16.0071    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   18.4792  -12.8252    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   16.8637  -12.1160    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   22.1353  -13.0395    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   22.1672  -13.9240    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   22.1809  -14.6581    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  5  1  1  0  0  0  0
  5  6  1  0  0  0  0
  7  3  1  0  0  0  0
  8  2  2  0  0  0  0
  9  3  1  0  0  0  0
 10  7  1  0  0  0  0
  7 11  1  0  0  0  0
  1 12  1  1  0  0  0
 13 14  1  0  0  0  0
 14  5  1  0  0  0  0
 15  4  2  0  0  0  0
 16 12  1  0  0  0  0
 17  8  1  0  0  0  0
 18 15  1  0  0  0  0
 14 19  1  6  0  0  0
 20 15  1  0  0  0  0
 21 11  1  0  0  0  0
 22 20  1  0  0  0  0
  3 23  1  1  0  0  0
  6  4  1  0  0  0  0
 16 11  1  0  0  0  0
  9 13  2  0  0  0  0
 10  8  1  0  0  0  0
 17 18  2  0  0  0  0
  2  1  1  0  0  0  0
  3  1  1  0  0  0  0
  4  2  1  0  0  0  0
 25 24  2  0  0  0  0
 26 24  1  0  0  0  0
 27 24  1  0  0  0  0
 28 24  2  0  0  0  0
  7 29  1  6  0  0  0
  5 30  1  1  0  0  0
 35 31  1  0  0  0  0
 35 36  1  6  0  0  0
 33 37  1  0  0  0  0
 38 32  2  0  0  0  0
 33 39  1  0  0  0  0
 37 40  1  0  0  0  0
 37 41  1  0  0  0  0
 31 42  1  0  0  0  0
 44 43  1  0  0  0  0
 35 44  1  0  0  0  0
 45 34  2  0  0  0  0
 46 42  1  0  0  0  0
 47 38  1  0  0  0  0
 48 45  1  0  0  0  0
 44 49  1  6  0  0  0
 50 45  1  0  0  0  0
 51 41  1  0  0  0  0
 52 50  1  0  0  0  0
 33 53  1  1  0  0  0
 36 34  1  0  0  0  0
 46 41  1  0  0  0  0
 39 43  2  0  0  0  0
 40 38  1  0  0  0  0
 47 48  2  0  0  0  0
 31 32  1  6  0  0  0
 31 33  1  0  0  0  0
 34 32  1  0  0  0  0
 37 54  1  6  0  0  0
M  END'''
        m = Chem.MolFromMolBlock(mb)
        self.assertEqual(len(Chem.GetMolFrags(m)), 6)
        smb, exclude = standardizer.get_parent_molblock(mb)
        sm = Chem.MolFromMolBlock(smb)
        self.assertEqual(len(Chem.GetMolFrags(sm)), 1)

    def testBondWedgingSurvivingSaltStripping(self):
        ' The salt stripper should not re-wedge bonds in the mol block '
        mb = '''sultamicillin
 

 42 46  0  0  1  0            999 V2000
    3.3125   -0.6500    0.0000 S   0  0
    3.3792    0.2583    0.0000 N   0  0
   -0.2333   -0.0875    0.0000 N   0  0
    2.8000   -0.4667    0.0000 C   0  0
    3.6750   -0.2042    0.0000 C   0  0  1  0  0  0
    2.8042    0.1125    0.0000 C   0  0  2  0  0  0
   -0.2333   -0.6500    0.0000 C   0  0  1  0  0  0
   -0.8083   -0.0875    0.0000 C   0  0
   -0.8083   -0.6500    0.0000 C   0  0  1  0  0  0
    3.8417    0.5625    0.0000 C   0  0
    0.3042    0.0875    0.0000 C   0  0  2  0  0  0
    0.3042   -0.8292    0.0000 S   0  0
    4.1417    0.0875    0.0000 C   0  0
    0.6417   -0.3792    0.0000 C   0  0
   -1.5208   -0.8167    0.0000 N   0  0
   -1.9708   -0.2875    0.0000 C   0  0
    2.5167    0.8458    0.0000 C   0  0
    0.7792    0.4208    0.0000 C   0  0
    3.0542   -1.2000    0.0000 O   0  0
    3.7917   -1.0542    0.0000 O   0  0
   -2.6500   -0.5250    0.0000 C   0  0  1  0  0  0
   -1.0250    0.3958    0.0000 O   0  0
    3.9542    1.1208    0.0000 O   0  0
   -1.9708    0.3083    0.0000 O   0  0
    2.0542    1.2708    0.0000 O   0  0
    1.4167    0.4958    0.0000 O   0  0
    1.5042    1.1000    0.0000 C   0  0
    3.0167    1.2208    0.0000 O   0  0
    0.5417    0.8500    0.0000 O   0  0
    2.4292   -0.9917    0.0000 C   0  0
    2.1542   -0.4042    0.0000 C   0  0
   -3.3708   -0.7250    0.0000 C   0  0
   -2.5125   -1.1375    0.0000 N   0  0
    0.9917   -0.8750    0.0000 C   0  0
    1.1292   -0.2292    0.0000 C   0  0
   -3.8458   -1.0042    0.0000 C   0  0
   -3.3250   -0.1542    0.0000 C   0  0
   -4.3458   -0.6917    0.0000 C   0  0
   -3.8250    0.1458    0.0000 C   0  0
   -4.3333   -0.1292    0.0000 C   0  0
   -0.2333   -1.2292    0.0000 H   0  0
    4.3000   -0.6417    0.0000 H   0  0
  2  5  1  0
  3 11  1  0
  4  1  1  0
  5  1  1  0
  6  4  1  0
  7  3  1  0
  8  3  1  0
  9  8  1  0
 10 13  1  0
 11 18  1  1
 12 14  1  0
 13  5  1  0
 14 11  1  0
  9 15  1  6
 16 15  1  0
  6 17  1  6
 18 26  1  0
 19  1  2  0
 20  1  2  0
 21 16  1  6
 22  8  2  0
 23 10  2  0
 24 16  2  0
 25 17  1  0
 26 27  1  0
 27 25  1  0
 28 17  2  0
 29 18  2  0
 30  4  1  0
 31  4  1  0
 21 32  1  6
 33 21  1  0
 34 14  1  0
 35 14  1  0
 36 32  1  0
 37 32  2  0
 38 36  2  0
 39 37  1  0
 40 39  2  0
  7 41  1  1
  5 42  1  6
  6  2  1  0
  2 10  1  0
  7 12  1  0
  7  9  1  0
 38 40  1  0
M  END'''
        smb, exclude = standardizer.get_parent_molblock(mb)
        inchi1 = Chem.MolBlockToInchi(mb)
        inchi2 = Chem.MolBlockToInchi(smb)
        self.assertEqual(inchi1,inchi2)


    def testParseMolblock(self):
        ' Check parse_molbock function '
        mb = '''
  Mrv1810 07121910172D          

  4  3  0  0  0  0            999 V2000
   -2.5038    0.4060    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -2.5038    1.2310    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -3.2182   -0.0065    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   -1.7893   -0.0065    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  4  0  0  0
M  CHG  2   2  -1   3   1
M  END'''
        mtt = standardizer.parse_molblock(mb, useRDKitChemistry=True)
        self.assertIs(type(mtt), Chem.Mol)
        mtf = standardizer.parse_molblock(mb, useRDKitChemistry=False)
        self.assertIs(type(mtf), Chem.Mol)


        mb = ''' wrong valence
  Mrv1810 07121910172D          

  4  3  0  0  0  0            999 V2000
   -2.5038    0.4060    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -2.5038    1.2310    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -3.2182   -0.0065    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   -1.7893   -0.0065    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  4  0  0  0
M  CHG  9   2  -1   3   1
M  END'''

        mft = standardizer.parse_molblock(mb, useRDKitChemistry=True)
        self.assertIsNone(mft)
        mff = standardizer.parse_molblock(mb, useRDKitChemistry=False)
        self.assertIsNone(mff)


        mb = '''
  MJ200400                      

 18 18  0  0  0  0  0  0  0  0999 V2000
   -0.8221    0.5653   -0.1449 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8063   -0.5926   -0.1512 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2643    0.9661   -0.1530 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2643   -0.9661   -1.1041 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8221   -0.5652   -1.1122 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8063    0.5926   -1.1059 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2211    0.8386   -0.8132 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2083    0.8314    0.5338 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4073   -0.8659    0.5171 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5804   -0.8782   -0.1453 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.2531    1.7911   -0.1485 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6634    0.6928    0.5152 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2531   -1.7911   -1.1086 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6634   -0.6928   -1.7724 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2211   -0.8386   -0.4439 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2083   -0.8313   -1.7910 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4074    0.8660   -1.7743 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.5804    0.8782   -1.1118 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  7  1  0  0  0  0
  1  8  1  0  0  0  0
  2  4  1  0  0  0  0
  2  9  1  0  0  0  0
  2 10  1  0  0  0  0
  3  6  1  0  0  0  0
  3 11  1  0  0  0  0
  3 12  1  0  0  0  0
  4  5  1  0  0  0  0
  4 13  1  0  0  0  0
  4 14  1  0  0  0  0
  5  6  1  0  0  0  0
  5 15  1  0  0  0  0
  5 16  1  0  0  0  0
  6 17  1  0  0  0  0
  6 18  1  0  0  0  0
M  END
'''

        ms = standardizer.parse_molblock(mb, useRDKitChemistry=True)
        mns = standardizer.parse_molblock(mb, useRDKitChemistry=False)
        self.assertNotEqual(ms.GetNumAtoms(), mns.GetNumAtoms())
