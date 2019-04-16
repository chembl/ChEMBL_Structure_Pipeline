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

    def test_triple_bonds(self):
        ms = [x for x in Chem.SDMolSupplier('./test_data/odd_drawings.sdf')]
        self.assertEqual(len(ms), 3)
        for m in ms:
            cm = standardizer.cleanup_drawing_mol(m)
            conf = cm.GetConformer()
            matches = m.GetSubstructMatches(Chem.MolFromSmarts('[#6]C#*'))
            self.assertTrue(len(matches))
            for match in matches:
                self.assertAlmostEqual(rdMolTransforms.GetAngleRad(
                    conf, match[0], match[1], match[2]), math.pi, places=2)
