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
    'Empty structure': 6,
    'Ambiguous stereo': 2
}


class InchiChecker(CheckerBase):
    name = "checks for InChI warnings"
    explanation = "checks for InChI warnings"
    penalty = 0
    element_matcher = re.compile(r'[A-Z][a-z]?.*?\([0-9]')

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
                return ((6, "InChI: Empty structure"),)
            else:
                if warnings2.find('Error 190 (no InChI; Unknown element(s)') == 0:
                    warnings2 = 'InChI: Unknown element(s)'
                return ((7, warnings2),)
        res = []
        for warning in warnings1.split(';'):
            warning = warning.strip()
            if warning == '':
                continue
            matched = False
            for k in inchiWarnings:
                if warning.find(k) == 0:
                    res.append((inchiWarnings[k], 'InChI: '+k))
                    matched = True
            if not matched:
                if InchiChecker.element_matcher.match(warning) and \
                        res[-1][1] == 'InChI: Accepted unusual valence(s)':
                    # The "Accepted unusual valence(s)" warning can produce things like this:
                    #  "Accepted unusual valence(s): N(4); Cu(6)"
                    # which break our splitting on ';' rule, producing a "warning" that's just
                    # element: "Cu(6)" in the example above.
                    # Here we detect this situation and avoid adding that element
                    # to the list of warnings
                    continue
                elif warning.find('bond(s)') == 0 and \
                        res[-1][1] == 'InChI: Ambiguous stereo':
                    # "Ambiguous stereo" can appear as:
                    # "Ambiguous stereo: center(s)", "Ambiguous stereo bond(s)" or
                    # "Ambiguous stereo: center(s); bond(s)"
                    # Here we detect the last case and avoid adding "bond(s)" to the
                    #  result
                    continue
                else:
                    res.append((2, 'InChI: '+warning))
        return tuple(sorted(res, reverse=True))


class StereoChecker(CheckerBase):
    name = "checks for stereo disagreements"
    explanation = "checks for stereo disagreements"
    penalty = 0
    stereo_matcher = re.compile('/t([^/]*)')

    @staticmethod
    def get_stereo_counts(molb):
        # InChI: count specified stereocenters in the *last* stereolayer
        nInchi = 0
        inchi, w1, w2 = get_inchi(molb)
        layers = StereoChecker.stereo_matcher.findall(inchi)
        if layers:
            nSpec = 0
            nUnspec = 0
            for layer in layers[-1].split(';'):
                if layer:
                    nSpecLayer = 0
                    nUnspecLayer = 0
                    n = 1
                    for center in layer.split(','):
                        center = center.replace(";", "")
                        if '*' in center:
                            n = int(re.sub('(\d+)\*.*','\\1',center))
                        if center[-1] == '?':
                            nUnspecLayer += 1
                        else:
                            nSpecLayer += 1
                    nSpec += nSpecLayer * n
            nInchi = int(nSpec)

        m = Chem.MolFromMolBlock(molb, sanitize=False, removeHs=False)

        # Mol blocks: count atoms where a wedged or hashed bond starts
        molCounter = Counter()
        for b in m.GetBonds():
            if b.HasProp('_MolFileBondStereo'):
                p = b.GetUnsignedProp('_MolFileBondStereo')
                if p == 1 or p == 6:
                    molCounter[b.GetBeginAtomIdx()] += 1
        nMol = len([k for k, v in molCounter.items() if v > 0])

        # RDKit molecule: count atoms that have tetrahedral stereo
        m.UpdatePropertyCache(False)
        Chem.AssignStereochemistry(m, force=True, cleanIt=True)
        nRDKit = len([1 for x in m.GetAtoms() if x.GetChiralTag() in (
            Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)])

        return nInchi, nMol, nRDKit

    @staticmethod
    def check(molb):
        """ returns true if there is any warning """
        nInchi, nMol, nRDKit = StereoChecker.get_stereo_counts(molb)
        if nInchi != nMol or nInchi != nRDKit:
            return True
        else:
            return False

    @staticmethod
    def get_stereo_score(molb):
        nInchi, nMol, nRDKit = StereoChecker.get_stereo_counts(molb)
        if nInchi != nMol and nInchi != nRDKit and nRDKit != nMol:
            return (5, 'Mol/Inchi/RDKit stereo mismatch')
        elif nRDKit == nMol and nInchi != nRDKit:
            return (5, 'RDKit_Mol/InChI stereo mismatch')
        elif nMol == nInchi and nInchi != nRDKit:
            return (2, 'InChi_Mol/RDKit stereo mismatch')
        elif nInchi == nRDKit and nInchi != nMol:
            return (5, 'InChi_RDKit/Mol stereo mismatch')
        return ()


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
                if abs(conf.GetAtomPosition(i).z) >= 0.0001:
                    return True
        return False


class Has3DFlagSetMolChecker(MolChecker):
    name = "has_3d_flag_set"
    explanation = "molecule has the 3D flag set for a 2D conformer"
    penalty = 2

    @staticmethod
    def check(mol):
        if mol.GetNumConformers():
            conf = mol.GetConformer()
            if conf.Is3D():
                for i in range(mol.GetNumAtoms()):
                    if abs(conf.GetAtomPosition(i).z) >= 0.0001:
                        return False

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
    penalty = 5

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

class HasManyOverlappingAtomsMolChecker(MolChecker):
    name = "has_many_overlapping_atoms"
    explanation = "molecule has six (or more) atoms with exactly the same coordinates"
    penalty = 6

    def check(mol):
        nOverlapping=0
        if mol.GetNumConformers():
            conf = mol.GetConformer()
            ps = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
            for i in range(len(ps)):
                for j in range(i):
                    d = ps[i]-ps[j]
                    if d.Length() < 0.0001:
                        nOverlapping += 1
                        if nOverlapping >= 6:
                            return True
        return False



class ZeroCoordsMolChecker(MolChecker):
    name = "zero_coordinates"
    explanation = "all atoms have zero coordinates"
    penalty = 6

    def check(mol):
        if not mol.GetNumAtoms() or mol.GetNumAtoms() == 1:
            return False
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
    explanation = "molecule has a radical that is not found in the known list"
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
             Has3DMolChecker, Has3DFlagSetMolChecker, HasIllegalBondTypeMolChecker, HasIllegalBondStereoMolChecker,
             HasMultipleStereoBondsMolChecker, HasManyOverlappingAtomsMolChecker, HasOverlappingAtomsMolChecker,
             ZeroCoordsMolChecker, HasCrossedRingBondMolChecker, HasStereoBondInRingMolChecker,
             HasStereoBondToStereocenterMolChecker, DisallowedRadicalMolChecker ]


def check_molblock(mb):
    mol = Chem.MolFromMolBlock(mb, sanitize=False, removeHs=False)
    if mol is None:
        return ((7, "Illegal input"),)
    res = []
    many_overlap = False
    for checker in _checkers:
        if issubclass(checker, MolFileChecker):
            matched = checker.check(mb)
        elif issubclass(checker, MolChecker):
            if checker.__name__ == 'HasOverlappingAtomsMolChecker' and many_overlap:
                continue
            matched = checker.check(mol)
            if checker.__name__ == 'HasManyOverlappingAtomsMolChecker' and matched:
                many_overlap = True
        else:
            raise ValueError(checker)
        if matched:
            res.append((checker.penalty, checker.explanation))
    res.extend(InchiChecker.get_inchi_score(mb))
    tpl = StereoChecker.get_stereo_score(mb)
    if tpl:
        res.append(tpl)
    return tuple(sorted(res, reverse=True))