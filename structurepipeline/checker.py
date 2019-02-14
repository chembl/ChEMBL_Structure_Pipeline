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
import rdkit

rdkversion = rdkit.__version__.split(".")[:2]
if rdkversion < ["2019", "03"]:
    raise ValueError("need an RDKit version >= 2019.03.1")


def _get_molblock_inchi_and_warnings(mb):
    inchi, res, w1, w2, auxinfo = rdinchi.MolBlockToInchi(mb)
    return inchi, w1, w2


class CheckerBase(object):
    __slots__ = ["name", "explanation", "penalty"]


class MolChecker(CheckerBase):
    pass


class MolFileChecker(CheckerBase):
    __slots__ = ["name", "explanation", "penalty"]


# used as a cache
__inchiDict = {}


def get_inchi(molb):
    h = hash(molb)
    if h not in __inchiDict:
        # make sure the cache doesn't get huge:
        if(len(__inchiDict) > 10000):
            __inchiDict.clear()
        __inchiDict[h] = _get_molblock_inchi_and_warnings(molb)
    return __inchiDict[h]


inchiWarnings = {
    'Accepted unusual valence(s)': 6,
    'Empty structure': 6
}


class InchiChecker(CheckerBase):
    name = "checks for InChI warnings"
    explanation = "checks for InChI warnings"
    penalty = 0

    @staticmethod
    def check(molb):
        """ returns true if there is any warning """
        inchi, w1, w2 = get_inchi(molb)
        if (not inchi) or w1 or w2:
            return True
        return False

    @staticmethod
    def get_inchi_score(molb):
        inchi, warnings1, warnings2 = get_inchi(molb)
        if not inchi:
            # no InChI was generated. This indicates either an error (penalty of 7) or an empty structure (penalty of 6)
            if molb.find('\n  0  0  ') > 0 and \
                    warnings2.find('Empty structure') > 0:
                return ((6, warnings2),)
            else:
                return ((7, warnings2),)
        res = []
        for warning in warnings1.split(';'):
            warning = warning.strip()
            if warning == '':
                continue
            matched = False
            for k in inchiWarnings:
                if warning.find(k) == 0:
                    res.append((inchiWarnings[k], k))
                    matched = True
            if not matched:
                res.append((2, warning))
        return tuple(sorted(res, reverse=True))


class NumAtomsMolChecker(MolChecker):
    name = "num_atoms_equals_zero"
    explanation = "number of atoms less than 1"
    penalty = 6

    @staticmethod
    def check(mol):
        return mol.GetNumAtoms() < 1


class Has3DMolChecker(MolChecker):
    name = "has_3d_coordinates"
    explanation = "molecule has 3D coordinates"
    penalty = 6

    @staticmethod
    def check(mol):
        if mol.GetNumConformers():
            conf = mol.GetConformer()
            for i in range(mol.GetNumAtoms()):
                if abs(conf.GetAtomPosition(i).z) > 0.0001:
                    return True
        return False


class HasIllegalBondTypeMolChecker(MolChecker):
    name = "has_illegal_bond_type"
    explanation = "molecule has a bond with an illegal type"
    penalty = 5

    @staticmethod
    def check(mol):
        legalTypes = (Chem.BondType.SINGLE,
                      Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)
        for bond in mol.GetBonds():
            if bond.GetBondType() not in legalTypes:
                return True
        return False


class HasIllegalBondStereoMolChecker(MolChecker):
    name = "has_illegal_bond_stereo"
    explanation = "molecule has a bond with an illegal stereo flag"
    penalty = 5

    @staticmethod
    def check(mol):
        for bond in mol.GetBonds():
            if bond.HasProp("_MolFileBondStereo") and \
                    bond.GetUnsignedProp("_MolFileBondStereo") == 4:
                return True
        return False


class HasMultipleStereoBondsMolChecker(MolChecker):
    name = "has_multiple_stereo_bonds"
    explanation = "molecule has an atom with multiple stereo bonds"
    penalty = 2

    @staticmethod
    def check(mol):
        atomsSeen = [0]*mol.GetNumAtoms()
        for bond in mol.GetBonds():
            if bond.HasProp("_MolFileBondStereo") and \
                    bond.GetUnsignedProp("_MolFileBondStereo"):
                if atomsSeen[bond.GetBeginAtomIdx()]:
                    return True
                atomsSeen[bond.GetBeginAtomIdx()] = 1
        return False


class HasOverlappingAtomsMolChecker(MolChecker):
    name = "has_overlapping_atoms"
    explanation = "molecule has two (or more) atoms with exactly the same coordinates"
    penalty = 6

    def check(mol):
        if mol.GetNumConformers():
            conf = mol.GetConformer()
            ps = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            for i in range(len(ps)):
                for j in range(i):
                    d = ps[i]-ps[j]
                    if d.Length() < 0.0001:
                        return True
        return False


class ZeroCoordsMolChecker(MolChecker):
    name = "zero_coordinates"
    explanation = "all atoms have zero coordinates"
    penalty = 6

    def check(mol):
        if mol.GetNumConformers():
            origin = Geometry.Point3D(0, 0, 0)
            conf = mol.GetConformer()
            for i in range(mol.GetNumAtoms()):
                p = conf.GetAtomPosition(i)
                d = p - origin
                if d.Length() > 0.0001:
                    return False
        return True


class HasCrossedRingBondMolChecker(MolChecker):
    name = "has_crossed_ring_bond"
    explanation = "molecule has a crossed bond in a ring"
    penalty = 5

    @staticmethod
    def check(mol):
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE and \
                    bond.GetBondDir() == Chem.BondDir.EITHERDOUBLE and \
                    bond.IsInRing():
                return True
        return False


class HasStereoBondInRingMolChecker(MolChecker):
    name = "has_stereobond_in_ring"
    explanation = "molecule has a stereo bond in a ring"
    penalty = 2

    @staticmethod
    def check(mol):
        for bond in mol.GetBonds():
            if bond.HasProp("_MolFileBondStereo") and \
                    bond.GetUnsignedProp("_MolFileBondStereo") != 3 and \
                    bond.IsInRing():
                return True
        return False


class HasStereoBondToStereocenterMolChecker(MolChecker):
    name = "has_stereo_bond_to_stereocenter"
    explanation = "molecule has an stereo bond to a stereocenter"
    penalty = 2

    @staticmethod
    def check(mol):
        atomsSeen = [0]*mol.GetNumAtoms()
        sbonds = []
        for bond in mol.GetBonds():
            if bond.HasProp("_MolFileBondStereo") and \
                    bond.GetUnsignedProp("_MolFileBondStereo"):
                atomsSeen[bond.GetBeginAtomIdx()] = 1
                sbonds.append(bond)
        for bond in sbonds:
            if atomsSeen[bond.GetEndAtomIdx()]:
                return True
        return False


class DisallowedRadicalMolChecker(MolChecker):
    name = "disallowed_radical"
    explanation = "molecule has a radical that isn't found in the known list"
    penalty = 6

    @staticmethod
    def check(mol):
        for atom in mol.GetAtoms():
            nrad = atom.GetNumRadicalElectrons()
            if(nrad > 0):
                if(atom.GetAtomicNum() == 7):
                    nbrs = atom.GetNeighbors()
                    # nitric oxide
                    if nrad != 1 or \
                            len(nbrs) != 1 or \
                            nbrs[0].GetAtomicNum() != 8 or \
                            mol.GetBondBetweenAtoms(atom.GetIdx(), nbrs[0].GetIdx()).GetBondType() != Chem.BondType.DOUBLE:
                        return True
                elif(atom.GetAtomicNum() == 8):
                    nbrs = atom.GetNeighbors()
                    # Aminoxyl
                    if nrad != 1 or \
                            len(nbrs) != 1 or \
                            nbrs[0].GetAtomicNum() != 7 or \
                            mol.GetBondBetweenAtoms(atom.GetIdx(), nbrs[0].GetIdx()).GetBondType() != Chem.BondType.SINGLE:
                        return True
                else:
                    return True
        return False


class PolymerFileChecker(MolFileChecker):
    name = "polymer_molfile"
    explanation = "polymer information in mol file"
    penalty = 6
    _regex = re.compile(
        r'^M  STY.+(SRU)|(MON)|(COP)|(CRO)|(ANY)', flags=re.MULTILINE)

    @staticmethod
    def check(data):
        return PolymerFileChecker._regex.search(data) is not None


class V3000FileChecker(MolFileChecker):
    name = "V3000_molfile"
    explanation = "V3000 mol file"
    penalty = 6
    _regex = re.compile(r'^M  [vV]30', flags=re.MULTILINE)

    @staticmethod
    def check(data):
        return V3000FileChecker._regex.search(data) is not None


_checkers = [PolymerFileChecker, V3000FileChecker, NumAtomsMolChecker,
             Has3DMolChecker, HasIllegalBondTypeMolChecker, HasIllegalBondStereoMolChecker,
             HasMultipleStereoBondsMolChecker, HasOverlappingAtomsMolChecker, ZeroCoordsMolChecker,
             HasCrossedRingBondMolChecker, HasStereoBondInRingMolChecker,
             HasStereoBondToStereocenterMolChecker, DisallowedRadicalMolChecker, ]


def check_molblock(mb):
    mol = Chem.MolFromMolBlock(mb, sanitize=False, removeHs=False)
    if mol is None:
        return (7, ((7, "Illegal input"),))
    res = []
    score = 0
    for checker in _checkers:
        if issubclass(checker, MolFileChecker):
            matched = checker.check(mb)
        elif issubclass(checker, MolChecker):
            matched = checker.check(mol)
        else:
            raise ValueError(checker)
        if matched:
            score += checker.penalty
            res.append((checker.penalty, checker.explanation))
    for penalty, msg in InchiChecker.get_inchi_score(mb):
        print(">>>>>>", penalty, msg)
        score += penalty
        res.append((penalty, msg))
    return score, tuple(sorted(res, reverse=True))
