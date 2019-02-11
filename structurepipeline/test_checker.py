#
#  Copyright (c) 2019 Greg Landrum
#  All rights reserved.
#
#  This file is part of the ChEMBL_StructurePipeline project.
#  The contents are covered by the terms of the MIT license
#  which is included in the file LICENSE, found at the root
#  of the source tree.

from . import checker
import unittest
from rdkit import Chem


class TestCase(unittest.TestCase):

    def test_num_atoms(self):
        m0 = Chem.Mol()
        m1 = Chem.MolFromSmiles('CC')
        self.assertFalse(checker.NumAtomsMolChecker.check(m1))
        self.assertTrue(checker.NumAtomsMolChecker.check(m0))

    def test_polymer(self):
        polyBlock = """
  MJ171200

 76 80  0  0  0  0  0  0  0  0999 V2000
   -6.4802    2.6494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8927    3.3638    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7177    3.3638    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1302    2.6494    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.7177    1.9349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8927    1.9349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2426    1.9349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6552    2.6494    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -6.4802    1.2203    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6552    1.2203    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0051    1.2203    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4176    1.9349    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2426    0.5060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4176    0.5060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8927    0.5060    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4802   -0.2085    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6552   -0.2085    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.1302    1.2203    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -6.8927   -0.9230    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4802   -1.6374    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6552   -1.6374    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0051   -1.6374    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.2426   -2.3519    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4176   -2.3519    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0051   -3.0663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1801   -3.0663    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7676   -3.7808    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4176   -3.7808    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0051   -4.4953    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1801   -4.4953    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3243   -3.7791    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7368   -3.0648    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5619   -3.0648    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9744   -3.7791    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5619   -4.4936    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7368   -4.4936    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3243   -5.2082    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9744   -5.2082    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5619   -5.9226    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7368   -5.9226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.9744   -2.3503    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7994   -2.3503    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7994   -3.7791    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.9487   -5.1497    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3243   -6.6371    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3705   -2.2987    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9580   -1.5842    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1329   -1.5842    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2796   -2.2987    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1329   -3.0132    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9580   -3.0132    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1329   -0.1553    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2796   -0.8698    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3705   -3.7276    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2796   -3.7276    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.1329   -4.4420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9580   -4.4420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3705   -5.1566    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8575   -2.2792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4450   -1.5648    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6200   -1.5648    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2073   -2.2792    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6200   -2.9937    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4450   -2.9937    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3789   -3.8003    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3797   -5.2304    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8445   -5.2489    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.3318   -2.2987    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.2075   -3.7082    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8195   -3.7288    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.3448   -2.2597    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1706   -0.8729    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3860   -0.6179    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    3.2997   -0.0580    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1046   -2.2987    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9675   -0.6233    0.0000 Na  0  3  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  1  6  1  0  0  0  0
  9 10  1  0  0  0  0
  7 10  1  0  0  0  0
  8  1  1  0  0  0  0
  6  9  1  0  0  0  0
 11 12  1  0  0  0  0
 13 14  1  0  0  0  0
 11 14  2  0  0  0  0
 10 13  2  0  0  0  0
 12  7  2  0  0  0  0
 15 16  2  0  0  0  0
 16 17  1  0  0  0  0
 13 17  1  0  0  0  0
  5 18  1  0  0  0  0
 19 20  1  0  0  0  0
 20 21  1  0  0  0  0
 16 19  1  0  0  0  0
 23 24  1  0  0  0  0
 22 24  2  0  0  0  0
 21 23  1  0  0  0  0
 25 26  1  0  0  0  0
 25 24  1  1  0  0  0
 28 29  1  0  0  0  0
 29 30  1  0  0  0  0
 27 30  1  0  0  0  0
 25 28  1  0  0  0  0
 27 26  1  0  0  0  0
 31 32  1  0  0  0  0
 32 33  1  0  0  0  0
 33 34  1  0  0  0  0
 34 35  1  0  0  0  0
 35 36  1  0  0  0  0
 31 36  1  0  0  0  0
 39 40  2  0  0  0  0
 37 40  1  0  0  0  0
 35 38  1  1  0  0  0
 36 37  1  6  0  0  0
 41 42  1  0  0  0  0
 34 43  1  6  0  0  0
 33 41  1  1  0  0  0
 44 38  1  1  0  0  0
 40 45  1  0  0  0  0
 46 47  1  0  0  0  0
 47 48  1  0  0  0  0
 48 49  1  0  0  0  0
 49 50  1  0  0  0  0
 50 51  1  0  0  0  0
 46 51  1  0  0  0  0
 52 53  1  0  0  0  0
 48 53  1  1  0  0  0
 56 57  2  0  0  0  0
 54 57  1  0  0  0  0
 51 54  1  6  0  0  0
 50 55  1  1  0  0  0
 57 58  1  0  0  0  0
 59 60  1  0  0  0  0
 60 61  1  0  0  0  0
 61 62  1  0  0  0  0
 62 63  1  0  0  0  0
 63 64  1  0  0  0  0
 59 64  1  0  0  0  0
 28 65  1  6  0  0  0
 31 65  1  1  0  0  0
 29 66  1  1  0  0  0
 30 67  1  6  0  0  0
 27 55  1  1  0  0  0
 46 68  1  1  0  0  0
 62 68  1  6  0  0  0
 63 69  1  1  0  0  0
 64 70  1  6  0  0  0
 59 71  1  1  0  0  0
 61 72  1  1  0  0  0
 72 73  1  0  0  0  0
 72 74  2  0  0  0  0
 49 75  1  6  0  0  0
M  STY  3   1 SRU   2 SRU   3 SRU
M  SCN  3   1 HT   2 HT   3 HT
M  SAL   1 15  55  50  51  54  57  58  56  46  68  62  63  69  64  70  61
M  SAL   1 11  72  74  73  60  47  48  53  52  49  75  59
M  SMT   1 b
M  SBL   1  2  71  76
M  SAL   2 15  27  26  25  28  65  31  36  37  40  45  39  35  34  43  33
M  SAL   2 15  41  42  32  29  66  24  22  23  21  20  19  16  17  13  10
M  SAL   2 15   7  12  11   9   6   1   8   2   3   4   5  18  14  15  30
M  SAL   2  2  67  38
M  SMT   2 a
M  SBL   2  2  71  46
M  SAL   3 15  38  35  36  37  40  45  39  31  65  28  25  24  22  23  21
M  SAL   3 15  20  19  16  17  13  10   7  12  11   9   6   1   8   2   3
M  SAL   3 15   4   5  18  14  15  26  27  55  50  51  54  57  58  56  46
M  SAL   3 15  68  62  63  69  64  70  61  72  74  73  60  47  48  53  52
M  SAL   3 13  49  75  30  67  29  66  32  33  41  42  34  43  59
M  SMT   3 n
M  SBL   3  2  46  76
M  END"""
        dataBlock = """

     RDKit          2D

  6  5  0  3  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8971    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1962   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4952    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  1  0
M  STY  3   1 MUL   2 SUP   3 DAT
M  SST  1   1 BLO
M  SCN  1   1 HH
M  SDS EXP  1   1
M  SPL  1   1   3
M  SNC  1   1   7
M  SBT  1   1   1
M  SAL   1  3   1   2   3
M  SPA   1  3   1   2   3
M  SBL   1  3   1   2   3
M  SDI   1  4    1.0000    3.0000    5.0000    7.0000
M  SDI   1  4    2.0000    4.0000    6.0000    8.0000
M  SMT   1 n
M  SBV   1   3
M  SAP   1  1   1   1 XX
M  SCL   1 TEST CLASS
M  SAL   2  3   4   5   6
M  SPA   2  3   4   5   6
M  SBL   2  3   3   4   5
M  SMT   2 TEST LABEL
M  SBV   2   3    3.0000    4.0000
M  SAP   2  1   4   0 YY
M  SDT   3 SAMPLE FIELD NAME              TSAMPLE FIELD INFO   PQSAMPLE QUERY OP
M  SDD   3 SAMPLE FIELD DISP
M  SED   3 SAMPLE DATA FIELD 1
M  SED   3 SAMPLE DATA FIELD 2
M  SED   3 SAMPLE DATA FIELD 3
M  END
$$$$
"""
        self.failUnless(checker.PolymerFileChecker.check(polyBlock))
        self.failIf(checker.PolymerFileChecker.check(dataBlock))

    def test_v3000(self):
        v3000Block = """
  Mrv1810 02111910072D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.75 0 0 0
M  V30 2 C 0.5837 0.77 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
        """
        v2000Block = """"
  Mrv1810 02111910062D          

  2  1  0  0  0  0            999 V2000
   -0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3127    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        self.failUnless(checker.V3000FileChecker.check(v3000Block))
        self.failIf(checker.V3000FileChecker.check(v2000Block))

    def test_has3d(self):
        matchBlock = """"
  Mrv1810 0211191006          

  2  1  0  0  0  0            999 V2000
   -0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3127    0.4125    0.1000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        nomatchBlock = """
  Mrv1810 0211191006

  2  1  0  0  0  0            999 V2000
   -0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3127    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        self.failUnless(checker.Has3DMolChecker.check(
            Chem.MolFromMolBlock(matchBlock, sanitize=False, removeHs=False)))
        self.failIf(checker.Has3DMolChecker.check(
            Chem.MolFromMolBlock(nomatchBlock, sanitize=False, removeHs=False)))

    def test_bondtypes(self):
        baseBlock = """
  Mrv1810 02111911082D          

  3  3  0  0  0  0            999 V2000
   -0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3127    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3333   -0.3745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  ZZZ  0  0  0  0
  3  2  1  0  0  0  0
M  END

"""
        for i in range(1, 9):
            mb = baseBlock.replace('ZZZ', str(i))
            m = Chem.MolFromMolBlock(mb, sanitize=False, removeHs=False)
            if i < 4:
                self.failIf(checker.HasIllegalBondTypeMolChecker.check(m))
            else:
                self.failUnless(checker.HasIllegalBondTypeMolChecker.check(m))

    def test_overlappingCoords(self):
        matchBlock = """"
  Mrv1810 0211191006          

  3  2  0  0  0  0            999 V2000
   -0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3127    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  3  2  1  0  0  0  0
M  END
"""
        nomatchBlock = """
  Mrv1810 0211191006

  3  2  0  0  0  0            999 V2000
   -0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3127    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  3  2  1  0  0  0  0
M  END
"""
        self.failUnless(checker.HasOverlappingAtomsMolChecker.check(
            Chem.MolFromMolBlock(matchBlock, sanitize=False, removeHs=False)))
        self.failIf(checker.HasOverlappingAtomsMolChecker.check(
            Chem.MolFromMolBlock(nomatchBlock, sanitize=False, removeHs=False)))

    def test_zeroCoords(self):
        matchBlock = """"
  Mrv1810 0211191006          

  2  1  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        nomatchBlock = """
  Mrv1810 0211191006

  2  1  0  0  0  0            999 V2000
   -0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3127    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        self.failUnless(checker.ZeroCoordsMolChecker.check(
            Chem.MolFromMolBlock(matchBlock, sanitize=False, removeHs=False)))
        self.failIf(checker.ZeroCoordsMolChecker.check(
            Chem.MolFromMolBlock(nomatchBlock, sanitize=False, removeHs=False)))

    def test_crossedRingbond(self):
        matchBlock = """"
  Mrv1810 02111912422D          

 13 13  0  0  0  0            999 V2000
    4.9010   -4.1237    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9002   -3.3190    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7011   -4.3045    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6966   -3.1458    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8407   -5.0870    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8447   -2.3664    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.6629   -5.1548    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.6654   -2.3039    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0091   -4.3755    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.0052   -3.0872    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.8252   -4.4613    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.8281   -2.9994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2722   -3.7315    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  2  4  1  0  0  0  0
  3  5  1  0  0  0  0
  4  6  1  0  0  0  0
  5  7  1  0  0  0  0
  6  8  1  0  0  0  0
  7  9  1  0  0  0  0
  8 10  1  0  0  0  0
  9 11  1  0  0  0  0
 10 12  1  0  0  0  0
 11 13  1  0  0  0  0
  1  2  1  0  0  0  0
 12 13  2  3  0  0  0
M  END
"""
        nomatchBlocks = ("""
  Mrv1810 02111912422D

 13 13  0  0  0  0            999 V2000
    4.9010   -4.1237    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9002   -3.3190    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7011   -4.3045    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6966   -3.1458    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8407   -5.0870    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8447   -2.3664    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.6629   -5.1548    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.6654   -2.3039    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0091   -4.3755    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.0052   -3.0872    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.8252   -4.4613    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.8281   -2.9994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2722   -3.7315    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  2  4  1  0  0  0  0
  3  5  1  0  0  0  0
  4  6  1  0  0  0  0
  5  7  1  0  0  0  0
  6  8  1  0  0  0  0
  7  9  1  0  0  0  0
  8 10  1  0  0  0  0
  9 11  1  0  0  0  0
 10 12  1  0  0  0  0
 11 13  1  0  0  0  0
  1  2  1  0  0  0  0
 12 13  2  0  0  0  0
M  END
""", """
  Mrv1810 02111912422D          

  4  3  0  0  0  0            999 V2000
    7.5223    0.2009    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.2368    0.6134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.9513    0.2009    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.6657    0.6134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  3  4  1  0  0  0  0
  2  3  2  3  0  0  0
M  END
""")
        self.failUnless(checker.HasCrossedRingBondMolChecker.check(
            Chem.MolFromMolBlock(matchBlock, sanitize=False, removeHs=False)))
        for mb in nomatchBlocks:
            self.failIf(checker.HasCrossedRingBondMolChecker.check(
                Chem.MolFromMolBlock(mb, sanitize=False, removeHs=False)))

    def test_disallowedRadicals(self):
        matches = [Chem.MolFromSmiles(x) for x in ('C[CH2]', 'C[O]', 'C[N]O')]
        nomatches = [Chem.MolFromSmiles(x) for x in ('[N]=O', 'N[O]', 'CN[O]')]
        for m in matches:
            self.failUnless(checker.DisallowedRadicalMolChecker.check(m))
        for m in nomatches:
            self.failIf(checker.DisallowedRadicalMolChecker.check(m))

    def test_illegalbondstereo(self):
        matchBlock = """"
  Mrv1810 0211191006          

  2  1  0  0  0  0            999 V2000
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  4  0  0  0
M  END
"""
        nomatchBlock = """
  Mrv1810 0211191006

  2  1  0  0  0  0            999 V2000
   -0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3127    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
M  END
"""
        self.failUnless(checker.HasIllegalBondStereoMolChecker.check(
            Chem.MolFromMolBlock(matchBlock, sanitize=False, removeHs=False)))
        self.failIf(checker.HasIllegalBondStereoMolChecker.check(
            Chem.MolFromMolBlock(nomatchBlock, sanitize=False, removeHs=False)))

    def test_multiplebondstereo(self):
        matchBlock = """"
  Mrv1810 02111915582D          

  5  4  0  0  0  0            999 V2000
   -2.7624    1.9506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0479    2.3631    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4768    2.3631    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -2.8293    1.1702    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5661    1.8059    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  6  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  1  2  1  1  0  0  0
M  END
"""
        nomatchBlock = """
  Mrv1810 02111915582D          

  5  4  0  0  0  0            999 V2000
   -2.7624    1.9506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0479    2.3631    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4768    2.3631    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -2.8293    1.1702    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5661    1.8059    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  1  2  1  1  0  0  0
M  END
"""
        self.failUnless(checker.HasMultipleStereoBondsMolChecker.check(
            Chem.MolFromMolBlock(matchBlock, sanitize=False, removeHs=False)))
        self.failIf(checker.HasMultipleStereoBondsMolChecker.check(
            Chem.MolFromMolBlock(nomatchBlock, sanitize=False, removeHs=False)))

        # another example, stereobonds *to* stereocenters do not count for this one:
        nomatchBlock = """
  Mrv1810 02111915582D          

  5  4  0  0  0  0            999 V2000
   -2.7624    1.9506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0479    2.3631    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4768    2.3631    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -2.8293    1.1702    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5661    1.8059    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  3  1  1  6  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  1  2  1  1  0  0  0
M  END
"""
        self.failIf(checker.HasMultipleStereoBondsMolChecker.check(
            Chem.MolFromMolBlock(nomatchBlock, sanitize=False, removeHs=False)))

    def test_stereobondInRing(self):
        matchBlock = """"
  Mrv1810 02111916052D          

  7  7  0  0  0  0            999 V2000
    3.3705    1.7732    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    2.7031    1.2883    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9580    0.5037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7830    0.5037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0380    1.2883    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1027    2.5536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0380    2.2581    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  1  5  1  0  0  0  0
  1  6  1  0  0  0  0
  1  7  1  0  0  0  0
  1  2  1  1  0  0  0
M  END
"""
        nomatchBlock = """
  Mrv1810 02111916062D          

  7  7  0  0  0  0            999 V2000
    3.3705    1.7732    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    2.7031    1.2883    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9580    0.5037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7830    0.5037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0380    1.2883    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.1027    2.5536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0380    2.2581    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  1  5  1  0  0  0  0
  1  6  1  0  0  0  0
  1  7  1  1  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        self.failUnless(checker.HasStereoBondInRingMolChecker.check(
            Chem.MolFromMolBlock(matchBlock, sanitize=False, removeHs=False)))
        self.failIf(checker.HasStereoBondInRingMolChecker.check(
            Chem.MolFromMolBlock(nomatchBlock, sanitize=False, removeHs=False)))

    def test_stereobondToStereocenter(self):
        matchBlock = """"
  Mrv1810 02111916122D          

  8  7  0  0  0  0            999 V2000
   -0.5357   -1.2277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1788   -0.8152    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -1.2502   -0.8152    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.5357   -2.0527    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3395   -1.4170    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    0.8932   -1.2277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1788    0.0098    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.9379   -0.5589    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  2  7  1  0  0  0  0
  2  8  1  0  0  0  0
  1  2  1  1  0  0  0
  2  6  1  1  0  0  0
M  END
"""
        nomatchBlock = """
  Mrv1810 02111916122D          

  8  7  0  0  0  0            999 V2000
   -0.5357   -1.2277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1788   -0.8152    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -1.2502   -0.8152    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.5357   -2.0527    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3395   -1.4170    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    0.8932   -1.2277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1788    0.0098    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.9379   -0.5589    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  1  0  0  0
  2  7  1  0  0  0  0
  2  8  1  0  0  0  0
  1  2  1  0  0  0  0
  2  6  1  1  0  0  0
M  END
"""
        self.failUnless(checker.HasStereoBondToStereocenterMolChecker.check(
            Chem.MolFromMolBlock(matchBlock, sanitize=False, removeHs=False)))
        self.failIf(checker.HasStereoBondToStereocenterMolChecker.check(
            Chem.MolFromMolBlock(nomatchBlock, sanitize=False, removeHs=False)))


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
