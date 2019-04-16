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
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdMolTransforms
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


# derived from the MolVS set
_normalization_transforms = """
//	Name	SMIRKS
Nitro to N+(O-)=O	[N;X3:1](=[O:2])=[O:3]>>[*+1:1]([*-1:2])=[*:3]
Sulfoxide to -S+(O-)-	[!O:1][S+0;D3:2](=[O:3])[!O:4]>>[*:1][S+1:2]([O-:3])[*:4]
Alkaline oxide to ions	[Li,Na,K;+0:1]-[O+0:2]>>([*+1:1].[O-:2])
Bad amide tautomer1	[C:1]([OH1;D1:2])=[NH1:3]>>[C:1](=[OH0:2])-[NH2:3]
Bad amide tautomer2	[C:1]([OH1;D1:2])=[NH0:3]>>[C:1](=[OH0:2])-[NH1:3]
"""
_normalizer_params = rdMolStandardize.CleanupParameters()
_normalizer = rdMolStandardize.NormalizerFromData(
    _normalization_transforms, _normalizer_params)


def normalize_mol(m):
    """
    """
    return _normalizer.normalize(m)


def remove_hs_from_mol(m):
    """ removes any "non-chiral" Hs

    Chiral Hs are:
    - Hs with a wedged/dashed bond to them

    """
    SENTINEL = 100
    for atom in m.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetDegree() == 1 and not atom.GetIsotope():
            nbr = atom.GetNeighbors()[0]
            bnd = atom.GetBonds()[0]
            if bnd.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH) or \
                    (bnd.HasProp("_MolFileBondStereo") and bnd.GetUnsignedProp("_MolFileBondStereo") in (1, 6)):
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


def get_fragment_parent_mol(m):
    with open('../data/solvents.smi') as inf:
        solvents = inf.read()
    solvent_remover = rdMolStandardize.FragmentRemoverFromData(
        solvents, skip_if_all_match=True)
    with open('../data/salts.smi') as inf:
        salts = inf.read()
    salt_remover = rdMolStandardize.FragmentRemoverFromData(
        salts, skip_if_all_match=True)
    res = salt_remover.remove(solvent_remover.remove(m))
    return res


def get_isotope_parent_mol(m):
    m = Chem.Mol(m)
    for at in m.GetAtoms():
        if at.GetIsotope():
            at.SetIsotope(0)
    return remove_hs_from_mol(m)


def get_parent_mol(m):
    return get_isotope_parent_mol(get_fragment_parent_mol(m))


def standardize_mol(m):
    m = update_mol_valences(m)
    m = remove_sgroups_from_mol(m)
    m = kekulize_mol(m)
    m = remove_hs_from_mol(m)
    m = normalize_mol(m)
    m = uncharge_mol(m)
    m = flatten_tartrate_mol(m)
    m = cleanup_drawing_mol(m)

    return m


def standardize_molblock(ctab):
    m = Chem.MolFromMolBlock(ctab, sanitize=False, removeHs=False)
    return Chem.MolToMolBlock(standardize_mol(m))
