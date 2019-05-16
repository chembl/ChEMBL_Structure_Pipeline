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

rdkversion = rdkit.__version__.split(".")[:2]
if rdkversion < ["2019", "03"]:
    raise ValueError("need an RDKit version >= 2019.03.1")


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
Sulfoxide to -S+(O-)-	[!O:1][S+0;D3:2](=[O:3])[!O:4]>>[*:1][S+1:2]([O-:3])[*:4]
Trivalent S	[O:1]=[S;D2;+0:2]-[#6:3]>>[*:1]=[*+1:2]-[*:3]
// Note that the next one doesn't work propertly because repeated appplications
// don't carry the cations from the previous rounds through. This should be
// fixed by implementing single-molecule transformations, but that's a longer-term
// project
//Alkaline oxide to ions	[Li,Na,K;+0:1]-[O+0:2]>>([*+1:1].[O-:2])
Bad amide tautomer1	[C:1]([OH1;D1:2])=[NH1:3]>>[C:1](=[OH0:2])-[NH2:3]
Bad amide tautomer2	[C:1]([OH1;D1:2])=[NH0:3]>>[C:1](=[OH0:2])-[NH1:3]
Halogen with no neighbors	[F,Cl,Br,I;X0;+0:1]>>[*-1:1]
"""
_normalizer_params = rdMolStandardize.CleanupParameters()
_normalizer = rdMolStandardize.NormalizerFromData(
    _normalization_transforms, _normalizer_params)

_alkoxide_pattern = Chem.MolFromSmarts('[Li,Na,K;+0]-[O+0]')


def normalize_mol(m):
    """
    """
    if m.HasSubstructMatch(_alkoxide_pattern):
        m = Chem.RWMol(m)
        for match in m.GetSubstructMatches(_alkoxide_pattern):
            m.RemoveBond(match[0], match[1])
            m.GetAtomWithIdx(match[0]).SetFormalCharge(1)
            m.GetAtomWithIdx(match[1]).SetFormalCharge(-1)
    return _normalizer.normalize(m)


def remove_hs_from_mol(m):
    """ removes most Hs

    Hs that are preserved by the RDKit's Chem.RemoveHs() will not
    be removed.

    Additional exceptions:
    - Hs with a wedged/dashed bond to them
    - Hs bonded to atoms with tetrahedral stereochemistry set
    - Hs bonded to atoms that have three (or more) ring bonds
    - Hs bonded to atoms with a charge != +1 in non-default valence states
    - Hs bonded to atoms with charge +1 in a valence state that's more than one 
      above the default

    """
    # we need ring info, so be sure it's there (this won't do anything if the rings
    # have already been found)
    Chem.FastFindRings(m)
    if m.NeedsUpdatePropertyCache():
        m.UpdatePropertyCache(strict=False)
    SENTINEL = 100
    for atom in m.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetDegree() == 1 and not atom.GetIsotope():
            nbr = atom.GetNeighbors()[0]
            bnd = atom.GetBonds()[0]
            preserve = False
            if bnd.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH) or \
                    (bnd.HasProp("_MolFileBondStereo") and bnd.GetUnsignedProp("_MolFileBondStereo") in (1, 6)):
                preserve = True
            else:
                if nbr.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW):
                    preserve = True
                elif (nbr.GetFormalCharge() != 1 and nbr.GetExplicitValence() > Chem.GetPeriodicTable().GetDefaultValence(nbr.GetAtomicNum())) or \
                        (nbr.GetFormalCharge() == 1 and nbr.GetExplicitValence() > Chem.GetPeriodicTable().GetDefaultValence(nbr.GetAtomicNum())+1):
                    preserve = True
                else:
                    ringBonds = [b for b in nbr.GetBonds() if
                                 m.GetRingInfo().NumBondRings(b.GetIdx())]
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
    return uncharger.uncharge(m)


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
    if(abs(abs(angle)-math.pi) > 0.017):
        rdMolTransforms.SetAngleRad(
            conf, nbrs[0], at.GetIdx(), nbrs[1], math.pi)


def _cleanup_triple_bonds(m):
    conf = m.GetConformer()
    if conf.Is3D():
        raise ValueError("can only operate on 2D conformers")
    for bond in m.GetBonds():
        if bond.GetBondType() == Chem.BondType.TRIPLE and m.GetRingInfo().NumBondRings(bond.GetIdx()) == 0:
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
        if(abs(abs(angle)-math.pi) > 0.017):
            rdMolTransforms.SetAngleRad(
                conf, match[0], match[1], match[2], math.pi)


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


_data_dir = os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "..", "data")
_solvents_file = os.path.join(_data_dir, "solvents.smi")
_salts_file = os.path.join(_data_dir, "salts.smi")


def get_fragment_parent_mol(m):
    basepath = os.path.dirname(os.path.abspath(__file__))
    with open(_solvents_file) as inf:
        solvents = inf.read()
    solvent_remover = rdMolStandardize.FragmentRemoverFromData(
        solvents, skip_if_all_match=True)
    with open(_salts_file) as inf:
        salts = inf.read()
    salt_remover = rdMolStandardize.FragmentRemoverFromData(
        salts, skip_if_all_match=True)
    # we need an aromatic representation for the salt/solvent removal to work, but
    # we don't want to lose the kekule form we came in with... so this is a bit messy
    m.UpdatePropertyCache(strict=False)
    nm = Chem.Mol(m)
    Chem.SetAromaticity(nm)
    tmp = salt_remover.remove(solvent_remover.remove(nm))
    if tmp.GetNumAtoms() != nm.GetNumAtoms():
        keep = set(nm.GetSubstructMatch(tmp))
        remove = set(range(m.GetNumAtoms())).difference(keep)
        res = Chem.RWMol(m)
        for idx in sorted(remove, reverse=True):
            res.RemoveAtom(idx)
        res = Chem.Mol(res)
    else:
        res = Chem.Mol(m)
    return res


def get_isotope_parent_mol(m):
    m = Chem.Mol(m)
    for at in m.GetAtoms():
        if at.GetIsotope():
            at.SetIsotope(0)
    return remove_hs_from_mol(m)


def get_parent_mol(m, neutralize=True):
    res = get_isotope_parent_mol(get_fragment_parent_mol(m))
    if neutralize:
        res = uncharge_mol(res)
    return res


def get_parent_molblock(ctab, neutralize=True):
    m = Chem.MolFromMolBlock(ctab, sanitize=False, removeHs=False)
    parent = get_parent_mol(m, neutralize=neutralize)
    return Chem.MolToMolBlock(parent)


def standardize_mol(m):
    exclude = exclude_flag(m, includeRDKitSanitization=False)
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
        if b.HasProp("_MolFileBondStereo"):
            val = b.GetProp("_MolFileBondStereo")
            if val == '1':
                b.SetBondDir(Chem.BondDir.BEGINWEDGE)
            elif val == '6':
                b.SetBondDir(Chem.BondDir.BEGINDASH)


def standardize_molblock(ctab):
    m = Chem.MolFromMolBlock(ctab, sanitize=False, removeHs=False)
    # the RDKit has, by default, removed bond wedging information from the molecule
    # put that back in:
    reapply_molblock_wedging(m)
    return Chem.MolToMolBlock(standardize_mol(m))
