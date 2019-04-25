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
import math
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdDepictor
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.info")


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
        data = (('[NH3+]CC[O-]', 'NCCO'),
                ('[NH3+]CCC[O-].[Na+]', 'NCCC[O-].[Na+]'),
                ('[Na+].[NH3+]CCC[O-]', 'NCCC[O-].[Na+]'),
                ('[NH3+]CCO', 'NCCO'),
                ('NCC[O-]', 'NCCO'),
                ('[Cl-].[NH3+]CCC[O-]', 'Cl.NCCCO'),
                ('[N+](C)(C)(C)CCC[O-]', '[N+](C)(C)(C)CCC[O-]'),
                ('[NH3+]CC([O-])C[O-]', 'NCC(O)CO'),
                # ('[NH3+]CC([O-])C[O-].[Na+]','NCC(O)C[O-].[Na+]'),
                ('[NH3+]CCC[O-].[NH+](C)(C)C', 'CN(C)C.NCCCO')
                )
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
        self.assertEqual(Chem.MolToSmiles(
            nm), "[K+].[K+].[K+].[K+].[K+].[K+].[O-]C1C([O-])C([O-])C([O-])C([O-])C1[O-]")
        nm = standardizer.standardize_mol(m)
        self.assertEqual(Chem.MolToSmiles(
            nm), "[K+].[K+].[K+].[K+].[K+].[K+].[O-]C1C([O-])C([O-])C([O-])C([O-])C1[O-]")

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
        tests = [('CN(C)(C)C', 'C[N+](C)(C)C'),
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
        self.assertEqual(sorted(standardizer._getAtomsToOtherSide(
            m.GetAtomWithIdx(3), m.GetBondBetweenAtoms(3, 4))), [0, 1, 2])
        self.assertEqual(sorted(standardizer._getAtomsToOtherSide(
            m.GetAtomWithIdx(4), m.GetBondBetweenAtoms(3, 4))), [5])

        m = Chem.MolFromSmiles('CC#N')
        self.assertEqual(sorted(standardizer._getAtomsToOtherSide(
            m.GetAtomWithIdx(1), m.GetBondBetweenAtoms(1, 2))), [0])
        self.assertEqual(sorted(standardizer._getAtomsToOtherSide(
            m.GetAtomWithIdx(2), m.GetBondBetweenAtoms(1, 2))), [])

    def test_triple_bonds_and_allenes(self):
        ms = [x for x in Chem.SDMolSupplier('./test_data/odd_drawings.sdf')]
        self.assertEqual(len(ms), 4)
        for m in ms:
            cm = standardizer.cleanup_drawing_mol(m)
            conf = cm.GetConformer()
            tmatches = m.GetSubstructMatches(Chem.MolFromSmarts('[#6]C#*'))
            amatches = m.GetSubstructMatches(Chem.MolFromSmarts('*C=C=C*'))
            self.assertTrue(len(tmatches) or len(amatches))
            for match in tmatches:
                self.assertAlmostEqual(rdMolTransforms.GetAngleRad(
                    conf, match[0], match[1], match[2]), math.pi, places=2)
            for match in amatches:
                self.assertAlmostEqual(rdMolTransforms.GetAngleRad(
                    conf, match[1], match[2], match[3]), math.pi, places=2)

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
        tests = [('c1cccnc1C(=O)O.[Na]', 'c1cccnc1C(=O)O'),
                 ('c1cccnc1C(=O)[O-].[Na+]', 'c1cccnc1C(=O)O'),
                 ('[Na].[Cl]', '[Na].[Cl]'), ('[Na+].[Cl-]', '[Na+].[Cl-]'),
                 ('c1cccnc1[NH3+].O=C(O)C(O)C(O)C(=O)O',
                  'c1cccnc1N'),
                 ('c1cccnc1[NH3+].O=C([O-])[C@H](O)[C@H](O)C(=O)[O-]',
                  'c1cccnc1N'),
                 ('c1cccnc1.ClCCl', 'c1cccnc1'),
                 ('c1cccnc1.ClCCl.[Na+].[Cl-].O', 'c1cccnc1'),
                 ('O=C([O-])C(O)C(O)C(=O)[O-]', 'O=C(O)C(O)C(O)C(=O)O'),
                 ('O=C([O-])C(O)C(O)C(=O)[O-].[Na+].[Na+]',
                  'O=C([O-])C(O)C(O)C(=O)[O-].[Na+].[Na+]'),
                 ]
        for smi, expected in tests:
            m = Chem.MolFromSmiles(smi)
            ssmi = Chem.MolToSmiles(standardizer.get_parent_mol(m))
            esmi = Chem.CanonSmiles(expected)
            self.assertEqual(ssmi, esmi)

    def test_fragment_parent2(self):
        tests = [('c1cccnc1C(=O)[O-].[Na+]', 'c1cccnc1C(=O)[O-]'),
                 ('[Na+].[Cl-]', '[Na+].[Cl-]'),
                 ('c1cccnc1[NH3+].O=C([O-])C(O)C(O)C(=O)[O-]',
                  'c1cccnc1[NH3+]'),
                 ('c1cccnc1[NH3+].O=C([O-])[C@H](O)[C@H](O)C(=O)[O-]',
                  'c1cccnc1[NH3+]'),
                 ('O=C([O-])C(O)C(O)C(=O)[O-]', 'O=C([O-])C(O)C(O)C(=O)[O-]'),
                 ('O=C([O-])C(O)C(O)C(=O)[O-].[Na+].[Na+]',
                  'O=C([O-])C(O)C(O)C(=O)[O-].[Na+].[Na+]'),
                 ]
        for smi, expected in tests:
            m = Chem.MolFromSmiles(smi)
            ssmi = Chem.MolToSmiles(
                standardizer.get_parent_mol(m, neutralize=False))
            esmi = Chem.CanonSmiles(expected)
            self.assertEqual(ssmi, esmi)
            # get_fragment_parent_mol doesn't do neutralization:
            ssmi = Chem.MolToSmiles(
                standardizer.get_fragment_parent_mol(m))

    def test_isotopes_parent1(self):
        tests = [('c1cc[13cH]nc1', 'c1cccnc1'),
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
        tests = [('c1cc[13cH]nc1.[Na]', 'c1cccnc1'),
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
            ssmi = Chem.MolToSmiles(standardizer.get_parent_mol(m))
            sp = Chem.SmilesParserParams()
            sp.removeHs = False
            em = Chem.MolFromSmiles(expected, sp)
            esmi = Chem.MolToSmiles(em)
            self.assertEqual(ssmi, esmi)

            # make sure we can also work from mol blocks:
            imb = Chem.MolToMolBlock(m)
            smb = standardizer.get_parent_molblock(imb)
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
        # This one we standardize
        self.assertFalse(mol is omol)
