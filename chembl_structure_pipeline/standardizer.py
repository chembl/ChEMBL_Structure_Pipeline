#
#  Copyright (c) 2019 Greg Landrum
#  All rights reserved.
#
#  This file is part of the ChEMBL_StructurePipeline project.
#  The contents are covered by the terms of the MIT license
#  which is included in the file LICENSE, found at the root
#  of the source tree.
import os
import re
from rdkit import Chem
from rdkit import Geometry
from rdkit.Chem import rdinchi
from collections import Counter
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdMolTransforms
from .exclude_flag import exclude_flag
import rdkit
import math
import sys

rdkversion = rdkit.__version__.split(".")
if rdkversion < ["2019", "09", "2"]:
    raise ValueError("need an RDKit version >= 2019.09.2")


def kekulize_mol(m):
    Chem.Kekulize(m)
    return m


def update_mol_valences(m):
    m = Chem.Mol(m)
    m.UpdatePropertyCache(strict=False)
    return m


# derived from the MolVS set, with ChEMBL-specific additions
_normalization_transforms = """
//	Name	SMIRKS
Nitro to N+(O-)=O	[N;X3:1](=[O:2])=[O:3]>>[*+1:1]([*-1:2])=[*:3]
Diazonium N	[*:1]-[N;X2:2]#[N;X1:3]>>[*:1]-[*+1:2]#[*:3]
Quaternary N	[N;X4;v4;+0:1]>>[*+1:1]
Trivalent O	[*:1]=[O;X2;v3;+0:2]-[#6:3]>>[*:1]=[*+1:2]-[*:3]
Sulfoxide to -S+(O-)	[!O:1][S+0;D3:2](=[O:3])[!O:4]>>[*:1][S+1:2]([O-:3])[*:4]
// this form addresses a pathological case that came up a few times in testing:
Sulfoxide to -S+(O-) 2	[!O:1][SH1+1;D3:2](=[O:3])[!O:4]>>[*:1][S+1:2]([O-:3])[*:4]
Trivalent S	[O:1]=[S;D2;+0:2]-[#6:3]>>[*:1]=[*+1:2]-[*:3]
// Note that the next one doesn't work propertly because repeated appplications
// don't carry the cations from the previous rounds through. This should be
// fixed by implementing single-molecule transformations, but that's a longer-term
// project
//Alkaline oxide to ions	[Li,Na,K;+0:1]-[O+0:2]>>([*+1:1].[O-:2])
Bad amide tautomer1	[C:1]([OH1;D1:2])=;!@[NH1:3]>>[C:1](=[OH0:2])-[NH2:3]
Bad amide tautomer2	[C:1]([OH1;D1:2])=;!@[NH0:3]>>[C:1](=[OH0:2])-[NH1:3]
Halogen with no neighbors	[F,Cl,Br,I;X0;+0:1]>>[*-1:1]
Odd pyridine/pyridazine oxide structure	[C,N;-;D2,D3:1]-[N+2;D3:2]-[O-;D1:3]>>[*-0:1]=[*+1:2]-[*-:3]
"""
_normalizer_params = rdMolStandardize.CleanupParameters()
_normalizer = rdMolStandardize.NormalizerFromData(_normalization_transforms,
                                                  _normalizer_params)

_alkoxide_pattern = Chem.MolFromSmarts('[Li,Na,K;+0]-[#7,#8;+0]')


def normalize_mol(m):
    """
    """
    Chem.FastFindRings(m)
    if m.HasSubstructMatch(_alkoxide_pattern):
        m = Chem.RWMol(m)
        for match in m.GetSubstructMatches(_alkoxide_pattern):
            m.RemoveBond(match[0], match[1])
            m.GetAtomWithIdx(match[0]).SetFormalCharge(1)
            m.GetAtomWithIdx(match[1]).SetFormalCharge(-1)
    res = _normalizer.normalize(m)
    return res


def remove_hs_from_mol(m):
    """ removes most Hs

    Hs that are preserved by the RDKit's Chem.RemoveHs() will not
    be removed.

    Additional exceptions:
    - Hs with a wedged/dashed bond to them
    - Hs bonded to atoms with tetrahedral stereochemistry set
    - Hs bonded to atoms that have three (or more) ring bonds that are not simply protonated
    - Hs bonded to atoms in a non-default valence state that are not simply protonated 


    For the above, the definition of "simply protonated" is an atom with charge = +1 and
    a valence that is one higher than the default.

    """
    # we need ring info, so be sure it's there (this won't do anything if the rings
    # have already been found)
    Chem.FastFindRings(m)
    if m.NeedsUpdatePropertyCache():
        m.UpdatePropertyCache(strict=False)
    SENTINEL = 100
    for atom in m.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetDegree(
        ) == 1 and not atom.GetIsotope():
            nbr = atom.GetNeighbors()[0]
            bnd = atom.GetBonds()[0]
            preserve = False
            if bnd.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH) or \
                    (bnd.HasProp("_MolFileBondStereo") and bnd.GetUnsignedProp("_MolFileBondStereo") in (1, 6)):
                preserve = True
            else:
                is_protonated = nbr.GetFormalCharge() == 1 and \
                    nbr.GetExplicitValence() == \
                    Chem.GetPeriodicTable().GetDefaultValence(nbr.GetAtomicNum())+1
                if nbr.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
                                          Chem.ChiralType.CHI_TETRAHEDRAL_CW):
                    preserve = True
                elif not is_protonated:
                    if nbr.GetExplicitValence() > Chem.GetPeriodicTable(
                    ).GetDefaultValence(nbr.GetAtomicNum()):
                        preserve = True
                    else:
                        ringBonds = [
                            b for b in nbr.GetBonds()
                            if m.GetRingInfo().NumBondRings(b.GetIdx())
                        ]
                        if len(ringBonds) >= 3:
                            preserve = True
            if preserve:
                # we're safe picking an arbitrary high value since you can't do this in a mol block:
                atom.SetIsotope(SENTINEL)

    res = Chem.RemoveHs(m, sanitize=False)
    for atom in res.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetIsotope() == SENTINEL:
            atom.SetIsotope(0)
    return res


def remove_sgroups_from_mol(m):
    """ removes all Sgroups
    """
    Chem.ClearMolSubstanceGroups(m)
    return m


def uncharge_mol(m):
    """

    >>> def uncharge_smiles(smi): return Chem.MolToSmiles(uncharge_mol(Chem.MolFromSmiles(smi)))
    >>> uncharge_smiles('[NH3+]CCC')
    'CCCN'
    >>> uncharge_smiles('[NH3+]CCC[O-]')
    'NCCCO'
    >>> uncharge_smiles('C[N+](C)(C)CCC[O-]')
    'C[N+](C)(C)CCC[O-]'
    >>> uncharge_smiles('CC[NH+](C)C.[Cl-]')
    'CCN(C)C.Cl'
    >>> uncharge_smiles('CC(=O)[O-]')
    'CC(=O)O'
    >>> uncharge_smiles('CC(=O)[O-].[Na+]')
    'CC(=O)[O-].[Na+]'
    >>> uncharge_smiles('[NH3+]CC(=O)[O-].[Na+]')
    'NCC(=O)[O-].[Na+]'
    >>> uncharge_smiles('CC(=O)[O-].C[NH+](C)C')
    'CC(=O)O.CN(C)C'

    Alcohols are protonated before acids:

    >>> uncharge_smiles('[O-]C([N+](C)C)CC(=O)[O-]')
    'C[N+](C)C(O)CC(=O)[O-]'

    And the neutralization is done in a canonical order, so atom ordering of the input
    structure isn't important:

    >>> uncharge_smiles('C[N+](C)(C)CC([O-])CC[O-]')
    'C[N+](C)(C)CC([O-])CCO'
    >>> uncharge_smiles('C[N+](C)(C)CC(CC[O-])[O-]')
    'C[N+](C)(C)CC([O-])CCO'

    """
    uncharger = rdMolStandardize.Uncharger(canonicalOrder=True)
    res = uncharger.uncharge(m)
    res.UpdatePropertyCache(strict=False)
    return res


def _getAtomsToOtherSide(startAt, bond):
    oAt = bond.GetOtherAtomIdx(startAt.GetIdx())
    res = []
    q = [x for x in startAt.GetNeighbors() if x.GetIdx() != oAt]
    while q:
        hd = q.pop(0)
        if hd.GetIdx() in res:
            continue
        res.append(hd.GetIdx())
        for nbr in hd.GetNeighbors():
            if nbr.GetIdx() == startAt.GetIdx():
                continue
            if nbr.GetIdx() == oAt:
                raise ValueError(f"cycle found {oAt} {res}")
            if nbr.GetIdx() not in res:
                q.append(nbr)
    return res


def _check_and_straighten_at_triple_bond(at, bond, conf):
    if at.GetDegree() != 2:
        raise ValueError("only works with degree 2")
    nbrs = [x.GetIdx() for x in at.GetNeighbors()]
    angle = rdMolTransforms.GetAngleRad(conf, nbrs[0], at.GetIdx(), nbrs[1])
    # are we off by more than a degree?
    if (abs(abs(angle) - math.pi) > 0.017):
        rdMolTransforms.SetAngleRad(conf, nbrs[0], at.GetIdx(), nbrs[1],
                                    math.pi)


def _cleanup_triple_bonds(m):
    conf = m.GetConformer()
    if conf.Is3D():
        raise ValueError("can only operate on 2D conformers")
    for bond in m.GetBonds():
        if bond.GetBondType() == Chem.BondType.TRIPLE and m.GetRingInfo(
        ).NumBondRings(bond.GetIdx()) == 0:
            at = bond.GetBeginAtom()
            if at.GetDegree() == 2:
                _check_and_straighten_at_triple_bond(at, bond, conf)
            at = bond.GetEndAtom()
            if at.GetDegree() == 2:
                _check_and_straighten_at_triple_bond(at, bond, conf)


def _cleanup_allenes(m):
    conf = m.GetConformer()
    if conf.Is3D():
        raise ValueError("can only operate on 2D conformers")
    p = Chem.MolFromSmarts('*=[C;R0]=*')
    for match in m.GetSubstructMatches(p):
        angle = rdMolTransforms.GetAngleRad(conf, match[0], match[1], match[2])
        # are we off by more than a degree?
        if (abs(abs(angle) - math.pi) > 0.017):
            rdMolTransforms.SetAngleRad(conf, match[0], match[1], match[2],
                                        math.pi)


def cleanup_drawing_mol(m):
    m = Chem.Mol(m)
    if not m.GetNumConformers():
        # if we don't have a conformer, just return
        return m
    conf = m.GetConformer()
    if conf.Is3D():
        for i in range(m.GetNumAtoms()):
            if abs(conf.GetAtomPosition(i).z) >= 0.0001:
                raise ValueError(
                    "cleanup_drawing_mol() only works for 2D molecules")
        conf.Set3D(False)
    Chem.FastFindRings(m)
    _cleanup_triple_bonds(m)
    _cleanup_allenes(m)
    return m


def flatten_tartrate_mol(m):
    tartrate = Chem.MolFromSmarts('OC(=O)C(O)C(O)C(=O)O')
    matches = m.GetSubstructMatches(tartrate)
    if matches:
        m = Chem.Mol(m)
        m.GetAtomWithIdx(3).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        m.GetAtomWithIdx(5).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    return m


_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
_solvents_file = os.path.join(_data_dir, "solvents.smi")
_salts_file = os.path.join(_data_dir, "salts.smi")


def get_fragment_parent_mol(m,
                            check_exclusion=False,
                            neutralize=False,
                            verbose=False):
    basepath = os.path.dirname(os.path.abspath(__file__))
    with open(_solvents_file) as inf:
        solvents = []
        for l in inf:
            if not l or l[0] == '#':
                continue
            l = l.strip().split('\t')
            if len(l) != 2:
                continue
            solvents.append((l[0], Chem.MolFromSmarts(l[1])))

    # there are a number of special cases for the ChEMBL salt stripping, so we
    # can't use the salt remover that's built into the RDKit standardizer.
    frags = []
    inputFrags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    for frag in inputFrags:
        frag = Chem.RemoveHs(frag, sanitize=False)
        frag.UpdatePropertyCache(strict=False)
        Chem.SetAromaticity(frag)
        frags.append(frag)
    keep = [1] * len(frags)
    for nm, solv in solvents:
        for i, frag in enumerate(frags):
            if keep[i] and frag.GetNumAtoms() == solv.GetNumAtoms() \
                and frag.GetNumBonds() == solv.GetNumBonds() \
                and frag.HasSubstructMatch(solv):
                keep[i] = 0
                if (verbose): print(f'matched solvent {nm}')
        if not max(keep):
            break
    if not max(keep):
        # everything removed, we can just return the input molecule:
        if check_exclusion:
            exclude = exclude_flag(m, includeRDKitSanitization=False)
        else:
            exclude = False
        if neutralize:
            res = uncharge_mol(m)
        else:
            res = Chem.Mol(m)
        return res, exclude

    with open(_salts_file) as inf:
        salts = []
        for l in inf:
            if not l or l[0] == '#':
                continue
            l = l.strip().split('\t')
            if len(l) != 2:
                continue
            salts.append((l[0], Chem.MolFromSmarts(l[1])))

    keepFrags1 = []
    keepFrags2 = []
    for i, v in enumerate(keep):
        if v:
            keepFrags1.append(frags[i])
            keepFrags2.append(inputFrags[i])
    frags = keepFrags1
    inputFrags = keepFrags2
    keep = [1] * len(frags)

    for nm, salt in salts:
        for i, frag in enumerate(frags):
            if keep[i] and frag.GetNumAtoms() == salt.GetNumAtoms() \
                and frag.GetNumBonds() == salt.GetNumBonds() \
                and frag.HasSubstructMatch(salt):
                if (verbose): print(f'matched salt {nm}')
                keep[i] = 0
        if not max(keep):
            break

    if not max(keep):
        # everything removed, keep everything:
        keep = [1] * len(frags)

    keepFrags = []
    seenSmis = set()
    for i, v in enumerate(keep):
        if not v:
            continue
        frag = inputFrags[i]
        if neutralize:
            cfrag = uncharge_mol(frag)
        else:
            cfrag = Chem.Mol(frag)
        keepFrags.append(i)
        # make sure there are no extraneous H atoms in the fragment:
        cfrag = Chem.RemoveHs(cfrag, sanitize=False)
        # need aromaticity perception to get a reasonable SMILES, but don't
        # want to risk a full sanitization:
        cfrag.ClearComputedProps()
        cfrag.UpdatePropertyCache(False)
        Chem.SanitizeMol(cfrag,
                         sanitizeOps=Chem.SANITIZE_SYMMRINGS
                         | Chem.SANITIZE_FINDRADICALS
                         | Chem.SANITIZE_SETAROMATICITY
                         | Chem.SANITIZE_ADJUSTHS)

        seenSmis.add(Chem.MolToSmiles(cfrag))
    if len(seenSmis) == 1:
        # if we just have one fragment left, this is easy:
        # just copy the fragment
        res = inputFrags[keepFrags[0]]
    else:
        # otherwise we need to create a molecule from the remaining fragments
        res = inputFrags[keepFrags[0]]
        for idx in keepFrags[1:]:
            frag = inputFrags[idx]
            res = Chem.CombineMols(res, frag)

    if check_exclusion:
        exclude = exclude_flag(res, includeRDKitSanitization=False)
    else:
        exclude = False

    # if we still match the exclude flag after stripping salts, go
    # back to the parent species after solvent stripping. These are now
    # in the inputFrags list
    if exclude:
        res = inputFrags[0]
        for frag in inputFrags[1:]:
            res = Chem.CombineMols(res, frag)

    if neutralize:
        res = uncharge_mol(res)
    return res, exclude


def get_isotope_parent_mol(m):
    m = Chem.Mol(m)
    for at in m.GetAtoms():
        if at.GetIsotope():
            at.SetIsotope(0)
    return remove_hs_from_mol(m)


def get_parent_mol(m, neutralize=True, check_exclusion=True, verbose=False):
    ipar = get_isotope_parent_mol(m)
    res, exclude = get_fragment_parent_mol(ipar,
                                           neutralize=neutralize,
                                           check_exclusion=check_exclusion,
                                           verbose=verbose)
    return res, exclude


def get_parent_molblock(ctab,
                        neutralize=True,
                        check_exclusion=True,
                        verbose=False):
    m = Chem.MolFromMolBlock(ctab, sanitize=False, removeHs=False)
    reapply_molblock_wedging(m)
    parent, exclude = get_parent_mol(m,
                                     neutralize=neutralize,
                                     check_exclusion=check_exclusion,
                                     verbose=verbose)
    return Chem.MolToMolBlock(parent, kekulize=False), exclude


def standardize_mol(m, check_exclusion=True):
    if check_exclusion:
        exclude = exclude_flag(m, includeRDKitSanitization=False)
    else:
        exclude = False
    if not exclude:
        m = update_mol_valences(m)
        m = remove_sgroups_from_mol(m)
        m = kekulize_mol(m)
        m = remove_hs_from_mol(m)
        m = normalize_mol(m)
        m = uncharge_mol(m)
        m = flatten_tartrate_mol(m)
        m = cleanup_drawing_mol(m)

    return m


def reapply_molblock_wedging(m):
    for b in m.GetBonds():
        # only do the wedgeing if the bond doesn't already have something there:
        if b.GetBondDir() == Chem.BondDir.NONE and b.HasProp(
                "_MolFileBondStereo"):
            val = b.GetProp("_MolFileBondStereo")
            if val == '1':
                b.SetBondDir(Chem.BondDir.BEGINWEDGE)
            elif val == '6':
                b.SetBondDir(Chem.BondDir.BEGINDASH)


def parse_molblock(ctab, useRDKitChemistry=False):
    if useRDKitChemistry:
        m = Chem.MolFromMolBlock(ctab, sanitize=True, removeHs=True)
    else:
        m = Chem.MolFromMolBlock(ctab, sanitize=False, removeHs=False)
        if not m:
            return None
        # the RDKit has, by default, removed bond wedging information from the molecule
        # put that back in:
        reapply_molblock_wedging(m)
        # Set the stereochemistry of double bonds
        # This block can be removed if github #X ends up being accepted and fixed
        anybonds = []
        for bond in m.GetBonds():
            if bond.GetStereo() == Chem.BondStereo.STEREOANY:
                anybonds.append(bond.GetIdx())
        Chem.SetBondStereoFromDirections(m)
        for bidx in anybonds:
            m.GetBondWithIdx(bidx).SetStereo(Chem.BondStereo.STEREOANY)
    return m


def standardize_molblock(ctab, check_exclusion=True):
    m = parse_molblock(ctab, useRDKitChemistry=False)
    if check_exclusion:
        if exclude_flag(m, includeRDKitSanitization=False):
            return ctab
    return Chem.MolToMolBlock(standardize_mol(m, check_exclusion=False))
