#
#  Copyright (c) 2019 Greg Landrum
#  All rights reserved.
#
#  This file is part of the ChEMBL_StructurePipeline project.
#  The contents are covered by the terms of the MIT license
#  which is included in the file LICENSE, found at the root
#  of the source tree.
from .. import checker
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
M  SCN  3   1 HT    2 HT    3 HT
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
        self.failUnlessEqual(checker.check_molblock(polyBlock),
                             ((6, 'polymer information in mol file'),
                              (5, 'InChi_RDKit/Mol stereo mismatch'),
                              (2, 'InChI: Proton(s) added/removed'),
                              (2, 'InChI: Ignore polymer data')))
        self.failUnlessEqual(checker.check_molblock(dataBlock),
                             ((2, 'InChI: Ignore polymer data'),))

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
        v2000Block = """
  Mrv1810 02111910062D

  2  1  0  0  0  0            999 V2000
   -0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3127    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        self.failUnless(checker.V3000FileChecker.check(v3000Block))
        self.failIf(checker.V3000FileChecker.check(v2000Block))
        self.failUnlessEqual(checker.check_molblock(v3000Block),
                             ((6, 'V3000 mol file'),))
        self.failUnlessEqual(checker.check_molblock(v2000Block),
                             ())

    def test_has3d(self):
        matchBlock = """
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
        self.failUnlessEqual(checker.check_molblock(matchBlock),
                             ((6, 'molecule has 3D coordinates'),))
        self.failUnlessEqual(checker.check_molblock(nomatchBlock),
                             ())

        matchBlock2 = """
  Mrv1810 02111910063D

  2  1  0  0  0  0            999 V2000
   -0.4018    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3127    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        self.failUnlessEqual(checker.check_molblock(matchBlock2),
                             ((2, 'molecule has the 3D flag set for a 2D conformer'),))

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
                self.failUnlessEqual(checker.check_molblock(mb),
                                     ())
        else:
            self.failUnless(checker.HasIllegalBondTypeMolChecker.check(m))
            self.failUnlessEqual(checker.check_molblock(mb),
                                 ((7, 'Error 108 (no InChI; Unrecognized bond type: 8) inp'),
                                  (5, 'molecule has a bond with an illegal type'),))

    def test_overlappingCoords(self):
        matchBlock = """
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
        self.failUnlessEqual(checker.check_molblock(matchBlock),
                             ((5, 'molecule has two (or more) atoms with exactly the same coordinates'),))
        self.failUnlessEqual(checker.check_molblock(nomatchBlock),
                             ())

    def test_zeroCoords(self):
        matchBlock = """
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
        self.failUnlessEqual(checker.check_molblock(matchBlock),
                             ((6, 'all atoms have zero coordinates'),
                             (5, 'molecule has two (or more) atoms with exactly the same coordinates')))
        self.failUnlessEqual(checker.check_molblock(nomatchBlock),
                             ())

    def test_crossedRingbond(self):
        matchBlock = """
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
""",)
        self.failUnless(checker.HasCrossedRingBondMolChecker.check(
            Chem.MolFromMolBlock(matchBlock, sanitize=False, removeHs=False)))
        self.failUnlessEqual(checker.check_molblock(matchBlock),
                             ((5, 'molecule has a crossed bond in a ring'),
                              (2, 'InChI: Omitted undefined stereo')))
        for mb in nomatchBlocks:
            self.failIf(checker.HasCrossedRingBondMolChecker.check(
                Chem.MolFromMolBlock(mb, sanitize=False, removeHs=False)))
            self.failUnlessEqual(checker.check_molblock(mb),
                                 ())

        nomatch2 = """
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
"""
        self.failIf(checker.HasCrossedRingBondMolChecker.check(
            Chem.MolFromMolBlock(nomatch2, sanitize=False, removeHs=False)))

        self.failUnlessEqual(checker.check_molblock(nomatch2),
                             ((2, 'InChI: Omitted undefined stereo'),))

    def test_disallowedRadicals(self):
        matches = [Chem.MolFromSmiles(x) for x in ('C[CH2]', 'C[O]', 'C[N]O')]
        nomatches = [Chem.MolFromSmiles(x) for x in ('[N]=O', 'N[O]', 'CN[O]')]
        for m in matches:
            self.failUnless(checker.DisallowedRadicalMolChecker.check(m))
            mb = Chem.MolToMolBlock(m)
            self.failUnlessEqual(checker.check_molblock(mb),
                                 ((6, "molecule has a radical that is not found in the known list"),))
        for m in nomatches:
            self.failIf(checker.DisallowedRadicalMolChecker.check(m))
            mb = Chem.MolToMolBlock(m)
            self.failUnlessEqual(checker.check_molblock(mb),
                                 ())

    def test_illegalbondstereo(self):
        matchBlock = """
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
        self.failUnlessEqual(checker.check_molblock(matchBlock),
                             ((5, 'molecule has a bond with an illegal stereo flag'),))
        self.failUnlessEqual(checker.check_molblock(nomatchBlock),
                             ((5, 'InChi_RDKit/Mol stereo mismatch'),))

    def test_multiplebondstereo(self):
        matchBlock = """
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
        self.failUnlessEqual(checker.check_molblock(matchBlock),
                             ((2, 'molecule has an atom with multiple stereo bonds'),))
        self.failUnlessEqual(checker.check_molblock(nomatchBlock),
                             ((2, 'InChI: Ambiguous stereo'),))

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
        self.failUnlessEqual(checker.check_molblock(nomatchBlock),
                             ((5, 'InChi_RDKit/Mol stereo mismatch'),
                              (2, 'molecule has an stereo bond to a stereocenter'),
                              (2, 'InChI: Ambiguous stereo')))

    def test_stereobondInRing(self):
        matchBlock = """
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
        self.failUnlessEqual(checker.check_molblock(matchBlock),
                             ((2, 'molecule has a stereo bond in a ring'),
                              (2, 'InChI: Ambiguous stereo')))
        self.failUnlessEqual(checker.check_molblock(nomatchBlock),
                             ())

    def test_stereobondToStereocenter(self):
        matchBlock = """
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
        self.failUnlessEqual(checker.check_molblock(matchBlock),
                             ((2, 'molecule has an stereo bond to a stereocenter'), (2, 'InChI: Ambiguous stereo')))
        self.failUnlessEqual(checker.check_molblock(nomatchBlock),
                             ())

    def test_inchiWarning(self):
        mb = """BDBM163075
     RDKit          2D

 27 30  0  0  0  0  0  0  0  0999 V2000
    1.6146   -5.5162    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9260   -4.0489    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0484    1.6535    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.2594    1.8470    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4379    3.0237    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1670    0.4398    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6489    0.4769    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3781    3.0608    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9796   -3.6396    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9460    3.1800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6752    0.5961    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1571    0.3205    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8863    2.9045    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0946   -2.6362    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4878   -3.7959    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5729    2.1226    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7839    1.3780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0647    1.9663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2758    1.5343    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0711   -1.5187    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6916    0.9089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6816   -0.1486    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4207   -1.6751    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1734    0.0078    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.1997    1.0652    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.3020   -0.4613    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.8110   -3.0456    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2 27  1  0
  3  5  2  0
  3  6  1  0
  4  7  2  0
  4  8  1  0
  5 10  1  0
  6 11  2  0
  7 12  1  0
  8 13  2  0
  9 14  2  0
  9 15  1  0
 10 18  2  0
 11 18  1  0
 12 19  2  0
 13 19  1  0
 14 20  1  0
 15 27  1  0
 16 18  1  0
 16 21  1  0
 17 19  1  0
 17 24  1  0
 20 22  2  0
 20 23  1  0
 21 25  2  0
 21 26  1  0
 22 24  1  0
 22 25  1  0
 23 26  2  0
 23 27  1  0
M  CHG  2  15  -1  27   1
M  END
"""
        self.failUnless(checker.InchiChecker.check(mb))
        r = checker.InchiChecker.get_inchi_score(mb)
        print(r)
        self.failUnlessEqual(
            r, ((6, 'InChI: Accepted unusual valence(s)'), (2, 'InChI: Charges were rearranged')))

    def test_inchiErrors(self):
        mb = """
  Mrv1810 02141913522D

  2  2  0  0  0  0            999 V2000
   -9.7074   -8.3118    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.7506   -9.1357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  1  1  0  0  0  0
M  END
"""
        self.failUnless(checker.InchiChecker.check(mb))
        r = checker.InchiChecker.get_inchi_score(mb)
        self.failUnlessEqual(
            r, ((7, 'Error 102 (no InChI; Multiple bonds between two atoms) inp'),))

        mb = """
  Mrv1810 02141913522D

  2  1  0  0  0  0            999 V2000
   -9.7074   -8.3118    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.7506   -9.1357    0.0000 Pz  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        self.failUnless(checker.InchiChecker.check(mb))
        r = checker.InchiChecker.get_inchi_score(mb)
        self.failUnlessEqual(
            r, ((7, 'InChI: Unknown element(s)'),))

        mb = """
  Mrv1810 02141914092D

  5  4  0  0  0  0            999 V2000
   -9.7074   -8.3118    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.7506   -9.1357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.5575   -9.3072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.5947   -9.8813    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.9256   -9.1357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  5  4  0  0  0  0
  2  4  4  0  0  0  0
  2  3  4  0  0  0  0
  1  2  4  0  0  0  0
M  END
"""
        self.failUnless(checker.InchiChecker.check(mb))
        r = checker.InchiChecker.get_inchi_score(mb)
        self.failUnlessEqual(
            r, ((7, 'Error 132 (no InChI; Atom has 1 or more than 3 aromatic bonds) inp'),))

        mb = """
  Mrv1810 02141913522D

  2  2  0  0  0  0            999 V2000
   -9.7074   -8.3118    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.7506   -9.1357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  6  1  0  0  0  0
M  END
"""
        self.failUnless(checker.InchiChecker.check(mb))
        r = checker.InchiChecker.get_inchi_score(mb)
        self.failUnlessEqual(
            r, ((7, 'Error 101 (no InChI; Bond to nonexistent atom) inp'),))

        mb = """
  Mrv1810 02141913522D

  0  0  0  0  0  0            999 V2000
M  END
"""
        self.failUnless(checker.InchiChecker.check(mb))
        r = checker.InchiChecker.get_inchi_score(mb)
        self.failUnlessEqual(
            r, ((6, 'InChI: Empty structure'),))
        r = checker.check_molblock(mb)
        self.failUnlessEqual(
            r, ((6, 'number of atoms less than 1'), (6, 'InChI: Empty structure'),))

        mb = """
  Mrv1810 02141914362D

 22 21  0  0  0  0            999 V2000
   -9.7074   -8.3118    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.7506   -9.1357    0.0000 Si  0  0  0  0  0  0  0  0  0  0  0  0
  -10.5575   -9.3072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.5947   -9.8813    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.9256   -9.1357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.3779   -8.5999    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0361   -9.5482    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.2698   -9.7768    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1094   -8.6165    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.5710   -9.0495    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.1672   -9.7191    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8796   -9.9505    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.0462   -8.3655    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.5352   -8.8808    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.4180   -9.6206    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.4150   -9.8894    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.9537   -8.9222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3012   -8.4438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.8796   -8.3208    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.2355   -8.4683    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -10.4856   -9.5102    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.9537   -9.3492    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  2  5  1  0  0  0  0
  2  6  1  0  0  0  0
  2  7  1  0  0  0  0
  2  8  1  0  0  0  0
  2  9  1  0  0  0  0
  2 10  1  0  0  0  0
  2 11  1  0  0  0  0
  2 12  1  0  0  0  0
  2 13  1  0  0  0  0
  2 14  1  0  0  0  0
  2 15  1  0  0  0  0
  2 16  1  0  0  0  0
  2 17  1  0  0  0  0
  2 18  1  0  0  0  0
  2 19  1  0  0  0  0
  2 20  1  0  0  0  0
  2 21  1  0  0  0  0
  2 22  1  0  0  0  0
M  END
"""
        self.failUnless(checker.InchiChecker.check(mb))
        r = checker.InchiChecker.get_inchi_score(mb)
        self.failUnlessEqual(
            r, ((7, "Error 104 (no InChI; Atom 'Si' has more than 20 bonds) inp"),))

    def test_illegalInput(self):
        mb = """

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
        self.failUnlessEqual(checker.check_molblock(mb),
                             ((7, 'Illegal input'),))

    def test_stereochecks(self):
        mb = """
  Mrv1810 02151908302D

 13 13  0  0  0  0            999 V2000
   -0.6027    0.3339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3171   -0.0786    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3171   -0.9036    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6027   -1.3161    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    0.1118   -0.9036    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1118   -0.0786    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -2.0316   -1.3161    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6027   -2.1411    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8263    0.3339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6027    1.1589    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    0.1118    1.5714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3171    1.5714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6920    1.9839    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  1  6  1  0  0  0  0
  3  7  1  0  0  0  0
  6  9  1  1  0  0  0
  4  8  1  1  0  0  0
  1 10  1  0  0  0  0
 10 11  1  0  0  0  0
 10 13  1  0  0  0  0
 10 12  1  1  0  0  0
M  ISO  1  12  13
M  END
"""
        self.assertEqual(checker.StereoChecker.get_stereo_counts(mb),
                         (3, 3, 3))
        self.assertFalse(checker.StereoChecker.check(mb))
        self.assertEqual(checker.StereoChecker.get_stereo_score(mb),
                         ())

        mb = """
  Mrv1719 01142015282D          

 53 56  0  0  0  0            999 V2000
    0.8250    0.0000    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    1.6500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0625   -0.7144    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4125   -0.7144    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125   -0.7144    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -0.8250   -1.4290    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6500   -1.4290    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0625   -2.1434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6500   -2.8579    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8250   -2.8579    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125   -2.1434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8250    0.0000    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -1.6500    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125    0.7144    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -0.8250    1.4290    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4125    0.7144    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    0.8250    1.4290    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.7295    0.0000    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    8.5545    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.9670   -0.7144    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.3170   -0.7144    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4920   -0.7144    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    6.0795   -1.4290    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2545   -1.4290    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8420   -2.1434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2545   -2.8579    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0795   -2.8579    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4920   -2.1434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0795    0.0000    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    5.2545    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4920    0.7144    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    6.0795    1.4290    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.3170    0.7144    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    7.7295    1.4290    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8250    7.2811    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    1.6500    7.2811    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0625    6.5666    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4125    6.5666    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125    6.5666    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -0.8250    5.8521    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6500    5.8521    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0625    5.1376    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6500    4.4232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8250    4.4232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125    5.1376    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8250    7.2811    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
   -1.6500    7.2811    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125    7.9955    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -0.8250    8.7100    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.4125    7.9955    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    0.8250    8.7100    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8875   -2.1434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0170   -2.1434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  1  0  0  0
  6  7  2  0  0  0  0
  7  8  1  0  0  0  0
  8  9  2  0  0  0  0
  9 10  1  0  0  0  0
 10 11  2  0  0  0  0
  5 12  1  0  0  0  0
 12 13  1  6  0  0  0
 12 14  1  0  0  0  0
 14 15  1  1  0  0  0
 14 16  1  0  0  0  0
 16 17  1  6  0  0  0
 18 19  1  1  0  0  0
 19 20  1  0  0  0  0
 18 21  1  0  0  0  0
 21 22  1  0  0  0  0
 22 23  1  1  0  0  0
 23 24  2  0  0  0  0
 24 25  1  0  0  0  0
 25 26  2  0  0  0  0
 26 27  1  0  0  0  0
 27 28  2  0  0  0  0
 22 29  1  0  0  0  0
 29 30  1  6  0  0  0
 29 31  1  0  0  0  0
 31 32  1  1  0  0  0
 31 33  1  0  0  0  0
 33 34  1  6  0  0  0
 35 36  1  1  0  0  0
 36 37  1  0  0  0  0
 35 38  1  0  0  0  0
 38 39  1  0  0  0  0
 39 40  1  1  0  0  0
 40 41  2  0  0  0  0
 41 42  1  0  0  0  0
 42 43  2  0  0  0  0
 43 44  1  0  0  0  0
 44 45  2  0  0  0  0
 39 46  1  0  0  0  0
 46 47  1  6  0  0  0
 46 48  1  0  0  0  0
 48 49  1  1  0  0  0
 48 50  1  0  0  0  0
 50 51  1  6  0  0  0
 16  1  1  0  0  0  0
 33 18  1  0  0  0  0
 50 35  1  0  0  0  0
 11  6  1  0  0  0  0
 28 23  1  0  0  0  0
 45 40  1  0  0  0  0
  8 52  1  0  0  0  0
 25 53  1  0  0  0  0
M  END
"""
        self.assertEqual(checker.StereoChecker.get_stereo_counts(mb),
                         (15, 15, 15))
        self.assertFalse(checker.StereoChecker.check(mb))
        self.assertEqual(checker.StereoChecker.get_stereo_score(mb),
                         ())

        mb = """

MJ200300                      
 21 18  0  0  1  0  0  0  0  0999 V2000
    1.4623   -0.2909    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1767   -0.7034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1767   -1.5284    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8912   -0.2909    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6057   -0.7034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6057   -1.5284    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3201   -0.2909    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0346   -0.7034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7491   -0.2909    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0346   -1.5284    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -0.5856   -0.3640    0.0000 Ca  0  2  0  0  0  0  0  0  0  0  0  0
    1.5521    2.3821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2665    1.9696    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2665    1.1446    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.9810    2.3821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6955    1.9696    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6955    1.1446    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.4099    2.3821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1244    1.9696    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8388    2.3821    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.1244    1.1446    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  2  3  1  6  0  0  0
  2  4  1  0  0  0  0
  5  4  1  0  0  0  0
  5  6  1  6  0  0  0
  5  7  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  2  0  0  0  0
  8 10  1  0  0  0  0
 13 12  1  0  0  0  0
 13 14  1  6  0  0  0
 13 15  1  0  0  0  0
 16 15  1  0  0  0  0
 16 17  1  6  0  0  0
 16 18  1  0  0  0  0
 18 19  1  0  0  0  0
 19 20  2  0  0  0  0
 19 21  1  0  0  0  0
M  END
"""

        self.assertEqual(checker.StereoChecker.get_stereo_counts(mb),
                          (4, 4, 4))
        self.assertFalse(checker.StereoChecker.check(mb))
        self.assertEqual(checker.StereoChecker.get_stereo_score(mb),
                          ())

        mb = """
  Mrv1810 02151908442D

  8  8  0  0  0  0            999 V2000
    5.2455    1.0482    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    4.5311    0.6357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5311   -0.1893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2455   -0.6018    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    5.9600   -0.1893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9600    0.6357    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2887    1.8721    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1593   -1.4223    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  1  6  1  0  0  0  0
  1  7  1  1  0  0  0
  4  8  1  1  0  0  0
M  END
"""
        self.assertEqual(checker.StereoChecker.get_stereo_counts(mb),
                         (2, 2, 2))
        self.assertFalse(checker.StereoChecker.check(mb))
        self.assertEqual(checker.StereoChecker.get_stereo_score(mb),
                         ())

        mb = """
  Mrv1810 02151908462D

  4  3  0  0  0  0            999 V2000
    2.2321    4.4196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0023    4.7153    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4117    4.5059    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9568    3.6420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END
"""
        self.assertEqual(checker.StereoChecker.get_stereo_counts(mb),
                         (0, 1, 0))
        self.assertTrue(checker.StereoChecker.check(mb))
        self.assertEqual(checker.StereoChecker.get_stereo_score(mb),
                         (5, 'InChi_RDKit/Mol stereo mismatch'))
        self.assertEqual(checker.check_molblock(mb),
                         ((5, 'InChi_RDKit/Mol stereo mismatch'),))

    def test_longermsg(self):
        mb = """CHEMBL3832892
  SciTegic03101709222D

 74 84  0  0  0  0            999 V2000
    0.2606    0.1503    0.0000 O   0  0
    1.3000    0.7500    0.0000 Cl  0  0
    2.3394    0.1503    0.0000 O   0  0
    1.3000    1.9500    0.0000 O   0  5
    1.3000   -0.4500    0.0000 O   0  0
    5.3605    0.1503    0.0000 O   0  0
    6.4000    0.7500    0.0000 Cl  0  0
    7.4394    0.1503    0.0000 O   0  0
    6.4000    1.9500    0.0000 O   0  5
    6.4000   -0.4500    0.0000 O   0  0
   14.5367   -8.3476    0.0000 O   0  0
   14.0883   -7.2345    0.0000 C   0  0
   12.6066   -7.0247    0.0000 N   0  0
   11.6871   -8.2111    0.0000 C   0  0
   10.4980   -8.0494    0.0000 C   0  0
   12.0293   -5.6344    0.0000 C   0  0
   10.8412   -5.4662    0.0000 O   0  0
   15.0108   -6.0564    0.0000 C   0  0
   16.4925   -6.2662    0.0000 C   0  0
   17.4330   -5.0907    0.0000 C   0  0
   16.8517   -3.7000    0.0000 C   0  0
   17.8559   -2.5847    0.0000 N   0  0
   17.5440   -1.1175    0.0000 N   0  0
   18.8418   -0.3693    0.0000 N   0  0
   19.9565   -1.3730    0.0000 C   0  0
   21.2555   -0.6230    0.0000 C   0  0
   20.9437    0.8442    0.0000 N   0  0
   22.1990    0.8690    0.0000 C   0  0
   21.8824    2.3352    0.0000 C   0  0
   20.7161    1.8113    0.0000 N   0  0
   21.2594    2.6508    0.0000 H   0  0
   20.7161    3.3113    0.0000 C   0  0
   19.4171    4.0612    0.0000 C   0  0
   18.1180    3.3112    0.0000 C   0  0
   18.1180    1.8112    0.0000 N   0  0
   18.2057    3.1738    0.0000 C   0  0
   19.6857    3.4175    0.0000 C   0  0
   20.7385    4.4838    0.0000 C   0  0
   22.0791    3.8110    0.0000 N   0  0
   23.4121    4.5005    0.0000 C   0  0
   24.6066    3.5959    0.0000 C   0  0
   25.9982    4.1702    0.0000 C   0  0
   26.1649    5.6664    0.0000 C   0  0
   27.5565    6.2408    0.0000 C   0  0
   27.7233    7.7370    0.0000 N   0  0
   29.0992    8.3366    0.0000 C   0  0
   29.2348    9.5290    0.0000 C   0  0
   28.5153    5.5192    0.0000 O   0  0
   24.9693    6.5663    0.0000 C   0  0
   23.5923    5.9810    0.0000 C   0  0
   22.3966    6.8809    0.0000 C   0  0
   22.5634    8.3771    0.0000 C   0  0
   23.9549    8.9514    0.0000 C   0  0
   25.1506    8.0515    0.0000 C   0  0
   26.5276    8.6369    0.0000 C   0  0
   26.6732    9.8280    0.0000 O   0  0
   21.8535    2.3281    0.0000 N   0  0
   20.3748    2.0851    0.0000 N   0  0
   19.4171    1.0612    0.0000 Cu  0  0
   17.0214    0.8177    0.0000 C   0  0
   17.6317   -0.5526    0.0000 C   0  0
   19.1235   -0.3955    0.0000 N   0  0
   18.6400   -1.2708    0.0000 H   0  0
   19.7519   -1.6748    0.0000 C   0  0
   21.2471   -1.5549    0.0000 C   0  0
   21.8909   -0.2001    0.0000 C   0  0
   19.3476   -2.7415    0.0000 C   0  0
   15.3740   -3.4906    0.0000 C   0  0
   14.8148   -2.1028    0.0000 C   0  0
   13.3331   -1.8930    0.0000 C   0  0
   12.3925   -3.0685    0.0000 C   0  0
   12.9698   -4.4589    0.0000 C   0  0
   14.4515   -4.6686    0.0000 C   0  0
   31.7688    0.7500    0.0000 O   0  0
  1  2  2  0
  2  3  2  0
  2  4  1  0
  2  5  2  0
  6  7  2  0
  7  8  2  0
  7  9  1  0
  7 10  2  0
 11 12  2  0
 12 13  1  0
 13 14  1  0
 14 15  1  0
 13 16  1  0
 16 17  2  0
 12 18  1  0
 18 19  2  3
 19 20  1  0
 20 21  2  3
 21 22  1  0
 22 23  1  0
 23 24  2  3
 24 25  1  0
 25 26  1  0
 26 27  1  0
 27 28  1  0
 28 29  1  0
 29 30  1  0
 30 31  1  0
 30 32  1  0
 32 33  1  0
 33 34  1  0
 34 35  1  0
 35 36  1  0
 36 37  1  0
 37 38  2  3
 38 39  1  0
 39 40  1  0
 40 41  2  3
 41 42  1  0
 42 43  2  3
 43 44  1  0
 44 45  1  0
 45 46  1  0
 46 47  1  0
 44 48  2  0
 43 49  1  0
 49 50  1  0
 50 40  1  0
 50 51  2  3
 51 52  1  0
 52 53  2  3
 49 54  2  3
 54 53  1  0
 54 55  1  0
 55 45  1  0
 55 56  2  0
 39 57  1  0
 57 58  2  3
 58 37  1  0
 58 59  1  0
 59 30  1  0
 59 27  1  0
 59 35  1  0
 59 24  1  0
 35 60  1  0
 60 61  1  0
 61 62  1  0
 62 59  1  0
 62 63  1  0
 62 64  1  0
 64 65  1  0
 65 66  1  0
 66 27  1  0
 25 67  2  3
 67 22  1  0
 21 68  1  0
 68 69  2  3
 69 70  1  0
 70 71  2  3
 71 72  1  0
 72 16  1  0
 72 73  2  3
 73 18  1  0
 73 68  1  0
M  CHG  2   4  -1   9  -1
M  END
"""
        self.assertEqual(checker.check_molblock(mb),
                         ((6, 'InChI: Accepted unusual valence(s)'),
                          (5, 'molecule has a crossed bond in a ring'),
                          (2, 'InChI: Proton(s) added/removed'),
                          (2, 'InChI: Metal was disconnected')))

        mb = """
  SciTegic11201309402D

 36 39  0  0  1  0            999 V2000
   -5.3728  -11.9974    0.0000 C   0  0
   -5.3728  -12.8224    0.0000 C   0  0
   -4.6590  -13.2288    0.0000 C   0  0  2  0  0  0
   -4.6590  -11.5788    0.0000 C   0  0
   -3.9495  -11.9645    0.0000 C   0  0
   -3.9488  -12.8224    0.0000 C   0  0  1  0  0  0
   -3.2355  -11.5836    0.0000 C   0  0  1  0  0  0
   -2.5205  -12.0033    0.0000 C   0  0  2  0  0  0
   -3.9526  -11.1724    0.0000 C   0  0
   -2.5220  -12.8280    0.0000 C   0  0  2  0  0  0
   -3.2440  -13.2301    0.0000 C   0  0  2  0  0  0
   -3.0808  -14.0406    0.0000 O   0  0
   -2.2603  -14.1399    0.0000 C   0  0
   -1.9160  -13.3883    0.0000 C   0  0
   -1.8585  -14.8569    0.0000 O   0  0
   -1.1066  -13.2252    0.0000 C   0  0
   -3.4803  -10.6975    0.0000 O   0  0
   -4.6703  -10.7607    0.0000 O   0  0
   -4.6737   -9.9357    0.0000 C   0  0
   -2.9511  -10.8481    0.0000 O   0  0
   -2.3191  -10.5772    0.0000 C   0  0
   -2.3236   -9.7522    0.0000 O   0  0
   -1.6308  -11.0304    0.0000 C   0  0
   -1.8142  -11.5828    0.0000 O   0  0
   -1.0179  -11.7963    0.0000 C   0  0
   -0.8045  -12.5925    0.0000 O   0  0
   -0.4310  -11.2134    0.0000 C   0  0  2  0  0  0
    0.3663  -11.4279    0.0000 C   0  0  1  0  0  0
    0.1478  -10.6265    0.0000 O   0  0
    0.9503  -12.0096    0.0000 C   0  0
   -0.8496  -10.4910    0.0000 C   0  0
   -3.9567  -13.6433    0.0000 O   0  0
   -5.0813  -13.9429    0.0000 C   0  0
   -4.0564  -13.8132    0.0000 H   0  0
   -1.6054  -12.4281    0.0000 H   0  0
   -3.4533  -13.6910    0.0000 H   0  0
 14 16  2  0
  8 10  1  0
  9 17  2  0
  2  3  1  0
  9 18  1  0
  3  6  1  0
 18 19  1  0
  5  9  1  0
  7 20  1  6
 10 11  1  0
 20 21  1  0
  5  4  2  0
 21 22  2  0
 21 23  1  0
  1  2  1  0
  8 24  1  1
  1  4  1  0
 24 25  1  0
  5  7  1  0
 25 26  2  0
 11 12  1  0
 25 27  1  0
 12 13  1  0
 27 28  1  0
 13 14  1  0
 28 29  1  0
 27 29  1  0
 10 14  1  0
 28 30  1  1
  6 11  1  0
 27 31  1  6
  3 32  1  0
 13 15  2  0
  6 32  1  0
  8  7  1  0
  3 33  1  1
 11 34  1  1
 10 35  1  6
  6 36  1  6
M  END
"""
        self.assertEqual(checker.check_molblock(mb),
                         ((2, 'InChI: Ambiguous stereo'),
                          ))

        mb = """713
  -OEChem-01301907122D

 24 27  0     0  0  0  0  0  0999 V2000
    3.8213    6.3167    0.0000 Tc  2  3  0  0  0  0  0  0  0  0  0  0
    2.3614    4.9700    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    3.8213    7.7167    0.0000 O   0  6  0  0  0  0  0  0  0  0  0  0
    5.6963    9.8099    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.5697    5.0663    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    7.5341    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.8337    2.1000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    7.6183    3.8739    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0886    7.5007    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.1964    4.9000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.3964    7.4900    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.7064    8.8200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3573    5.7663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    6.8341    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3573    6.8665    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9914    8.8200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1964    2.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    5.7663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4088    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9839    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.6212    2.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4088    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9839    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1964    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  2  0  0  0  0
  1  9  1  0  0  0  0
  1 10  1  0  0  0  0
  1 11  1  0  0  0  0
  2 18  1  0  0  0  0
  4 12  2  0  0  0  0
  5 13  2  0  0  0  0
  6 14  2  0  0  0  0
  7 21  1  0  0  0  0
  8 21  2  0  0  0  0
  9 12  1  0  0  0  0
  9 15  1  0  0  0  0
 10 13  1  0  0  0  0
 10 17  1  0  0  0  0
 11 14  1  0  0  0  0
 11 16  1  0  0  0  0
 12 16  1  0  0  0  0
 13 15  1  0  0  0  0
 14 18  1  0  0  0  0
 17 19  1  0  0  0  0
 17 20  2  0  0  0  0
 19 21  1  0  0  0  0
 19 22  2  0  0  0  0
 20 23  1  0  0  0  0
 22 24  1  0  0  0  0
 23 24  2  0  0  0  0
M  CHG  3   1   1   3  -2   7  -1
M  ISO  1   1  99
M  END"""
        self.assertEqual(checker.check_molblock(mb),
                         ((6, 'InChI: Accepted unusual valence(s)'),
                          (2, 'InChI: Proton(s) added/removed'),
                          (2, 'InChI: Metal was disconnected')))

        mb = """'237894132
  -OEChem-01301913592D

107111  0     1  0  0  0  0  0999 V2000
   -0.2805   -6.2347    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0763   -6.2347    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
   -8.7592   -4.5847    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
    4.2274   -6.2347    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
   11.2772   -6.2347    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
   13.5312   -6.2347    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
   12.7062   -8.7097    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
  -11.0103   -1.6252    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -11.7247   -0.3877    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -10.2958   -2.8627    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -11.0103   -1.6252    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -11.4860   -6.0076    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -13.6066   -2.7799    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -14.1786   -5.7543    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -14.8441   -2.0655    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -13.1682   -4.7439    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3935   -7.0597    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6474   -5.4097    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.4737   -5.8222    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3138   -6.9492    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1485   -5.4097    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.4024   -7.0597    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.9103   -7.0597    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9329   -9.1222    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.6563   -5.4097    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.4522   -7.0597    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   12.7062   -5.4097    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -12.4392   -2.4502    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  -14.0191   -3.4944    0.0000 N   0  0  3  0  0  0  0  0  0  0  0  0
  -13.7516   -4.1605    0.0000 N   0  0  3  0  0  0  0  0  0  0  0  0
   -0.1239   -4.1722    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  -10.9026   -4.9972    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8224   -6.2347    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.9735   -6.2347    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    9.0233   -6.2347    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    6.4813   -6.2347    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    7.9103   -3.7597    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    9.1798   -4.1722    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3619   -9.1222    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  -12.6527   -3.2471    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.3672   -3.6596    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.1094   -4.3644    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -12.0694   -3.8305    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -12.4392   -1.6252    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.7247   -2.8627    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -14.7088   -4.1744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -14.4220   -4.8212    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -12.2829   -4.6273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.9651   -4.9574    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -14.4316   -2.7799    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.0103   -2.4502    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.7247   -1.2127    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3935   -3.7597    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3935   -4.5847    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -11.6995   -5.2107    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6088   -3.5047    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -14.0191   -2.0655    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -13.3817   -5.5408    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6088   -4.8396    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1079   -4.9972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1079   -5.8222    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -2.0065   -3.2076    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4373   -2.6978    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3935   -6.2347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8350   -2.4007    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0504   -2.1457    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6474   -6.2347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3619   -6.6472    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
  -10.1882   -4.5847    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4340   -6.6472    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -9.4737   -4.9972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9013   -6.2347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6158   -5.8222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1485   -6.2347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.0447   -4.9972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3619   -7.4722    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6879   -5.8222    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -7.3303   -6.2347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.6158   -4.9972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.0447   -5.8222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3303   -4.5847    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4024   -6.2347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4340   -7.4722    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6474   -7.8847    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6879   -4.9972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1958   -5.8222    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
    7.9103   -6.2347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.7377   -5.8222    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -3.6474   -8.7097    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1958   -4.9972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6563   -6.2347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4522   -6.2347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.9103   -4.5847    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9419   -6.6472    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.7377   -4.9972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4024   -4.5847    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9735   -4.5847    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.9917   -6.6472    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
    8.6949   -4.8396    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4522   -4.5847    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.6949   -3.5047    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.7062   -6.2347    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.9917   -7.4722    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.1667   -4.9972    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.4522   -3.7597    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.7062   -7.8847    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   13.4206   -9.1222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1 64  1  0  0  0  0
  1 70  1  0  0  0  0
  2 68  1  0  0  0  0
  2 72  1  0  0  0  0
  3 71  1  0  0  0  0
  3 75  1  0  0  0  0
  4 82  1  0  0  0  0
  4 94  1  0  0  0  0
  5 92  1  0  0  0  0
 98  5  1  6  0  0  0
  6102  1  0  0  0  0
  7106  1  0  0  0  0
  7107  1  0  0  0  0
  8 51  2  0  0  0  0
  9 52  2  0  0  0  0
 10 51  1  0  0  0  0
 11 52  1  0  0  0  0
 12 55  2  0  0  0  0
 13 57  2  0  0  0  0
 14 58  2  0  0  0  0
 15 57  1  0  0  0  0
 16 58  1  0  0  0  0
 17 64  2  0  0  0  0
 18 67  2  0  0  0  0
 19 71  2  0  0  0  0
 20 72  2  0  0  0  0
 21 74  2  0  0  0  0
 22 82  2  0  0  0  0
 23 87  2  0  0  0  0
 24 89  2  0  0  0  0
 25 91  2  0  0  0  0
 26 92  2  0  0  0  0
 27102  2  0  0  0  0
 28 40  1  0  0  0  0
 28 44  1  0  0  0  0
 28 45  1  0  0  0  0
 29 41  1  0  0  0  0
 29 46  1  0  0  0  0
 29 50  1  0  0  0  0
 30 42  1  0  0  0  0
 30 47  1  0  0  0  0
 30 49  1  0  0  0  0
 31 56  1  0  0  0  0
 31 59  1  0  0  0  0
 32 55  1  0  0  0  0
 32 69  1  0  0  0  0
 61 33  1  1  0  0  0
 33 67  1  0  0  0  0
 34 74  1  0  0  0  0
 77 34  1  1  0  0  0
 35 87  1  0  0  0  0
 88 35  1  1  0  0  0
 86 36  1  1  0  0  0
 36 91  1  0  0  0  0
 37 93  1  0  0  0  0
 37101  1  0  0  0  0
 38 99  1  0  0  0  0
 38101  2  0  0  0  0
 39 89  1  0  0  0  0
 40 41  1  0  0  0  0
 40 42  1  0  0  0  0
 40 43  1  0  0  0  0
 43 48  1  0  0  0  0
 44 52  1  0  0  0  0
 45 51  1  0  0  0  0
 46 47  1  0  0  0  0
 48 55  1  0  0  0  0
 49 58  1  0  0  0  0
 50 57  1  0  0  0  0
 53 54  1  0  0  0  0
 53 56  2  0  0  0  0
 53 62  1  0  0  0  0
 54 59  2  0  0  0  0
 54 60  1  0  0  0  0
 56 63  1  0  0  0  0
 60 61  1  0  0  0  0
 61 64  1  0  0  0  0
 62 65  2  0  0  0  0
 63 66  2  0  0  0  0
 65 66  1  0  0  0  0
 68 67  1  6  0  0  0
 68 76  1  0  0  0  0
 69 71  1  0  0  0  0
 70 74  1  0  0  0  0
 70 83  1  6  0  0  0
 72 73  1  0  0  0  0
 73 78  2  0  0  0  0
 73 79  1  0  0  0  0
 75 80  2  0  0  0  0
 75 81  1  0  0  0  0
 76 84  1  0  0  0  0
 77 82  1  0  0  0  0
 77 85  1  0  0  0  0
 78 80  1  0  0  0  0
 79 81  2  0  0  0  0
 84 89  1  0  0  0  0
 85 96  1  0  0  0  0
 85 97  1  0  0  0  0
 86 87  1  0  0  0  0
 86 90  1  0  0  0  0
 88 92  1  0  0  0  0
 88 95  1  0  0  0  0
 90 93  1  0  0  0  0
 91 94  1  0  0  0  0
 93 99  2  0  0  0  0
 95100  1  0  0  0  0
 98102  1  0  0  0  0
 98103  1  0  0  0  0
100104  1  0  0  0  0
100105  1  0  0  0  0
103106  1  0  0  0  0
A    1
H
A    2
H
A    3
H
A    4
H
A    5
H
A    6
NH2.
M  END"""
        self.assertEqual(checker.check_molblock(mb),
                         ((6, 'InChI: Accepted unusual valence(s)'),
                         (5, 'molecule has two (or more) atoms with exactly the same coordinates')))

    def test_MultipleOVerlapping(self):
        mb = """
  Mrv1824 10091914432D          

 14  7  0  0  0  0            999 V2000
   -1.9420    1.3161    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6564    0.9036    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6564    0.0786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9420   -0.3339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2275    0.0786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2275    0.9036    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9420    2.1411    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9420    1.3161    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6564    0.9036    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6564    0.0786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9420   -0.3339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2275    0.0786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2275    0.9036    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9420    2.1411    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  1  6  2  0  0  0  0
  1  7  1  0  0  0  0
M  END"""
        self.assertEqual(checker.check_molblock(mb),
  ((6, 'molecule has six (or more) atoms with exactly the same coordinates'),))

