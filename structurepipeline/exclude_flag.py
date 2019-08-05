from rdkit import Chem
import os

_data_dir = os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "data")

# list of mols to keep no matter what
mols_to_keep = ['exclude_mols/ranitidine_to_keep.mol', 'exclude_mols/CO-ADD1.mol',
                'exclude_mols/CO-ADD2.mol', 'exclude_mols/CO-ADD3.mol', 'exclude_mols/CO-ADD4.mol']

# structures with vanadium that we want to keep
vanadium_to_keep = ['exclude_mols/v1.mol', 'exclude_mols/v2.mol',
                    'exclude_mols/v3.mol', 'exclude_mols/v4.mol']

# Zn not in list mimicking old PP protocol
METAL_LIST = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Ga', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
              'Cd', 'In', 'Sn', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Ac',
              'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Th', 'Pa', 'U', 'Np',
              'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Ge', 'Sb']


def exclude_flag(molfile, dataDir='', includeRDKitSanitization=True):
    """
    Rules to exclude structures. All exceptions are because of the old PP protocol.
    """
    if not dataDir:
        dataDir = _data_dir
    rdkit_fails = False
    exclude = False
    metallic = False
    boron_count = 0
    c_count = 0
    ok_count = 0

    if type(molfile) == str:
        mol = Chem.MolFromMolBlock(molfile, sanitize=False)
        if includeRDKitSanitization:
            try:
                Chem.SanitizeMol(mol)
            except:
                rdkit_fails = True
    else:
        mol = molfile
    # exclude all molecules that can't be sanitized by RDKit except the ones in the two weird sets

    # molecules that won't be excluded, mols_to_keep
    matches = [mol.HasSubstructMatch(Chem.MolFromMolFile(
        os.path.join(dataDir, x), sanitize=False)) for x in mols_to_keep]
    if True in matches:
        return exclude

    # structures with vanadium to keep
    matches = [mol.HasSubstructMatch(Chem.MolFromMolFile(
        os.path.join(dataDir, x), sanitize=False)) for x in vanadium_to_keep]
    # if any saved vanadium containing molecule matches with NotV, will be also excluded
    if True in matches and not mol.HasSubstructMatch(Chem.MolFromMolFile(os.path.join(dataDir, 'exclude_mols/NotV.mol'), sanitize=False)):
        return exclude

    for atom in mol.GetAtoms():
        a_type = atom.GetSymbol()
        if a_type in METAL_LIST:
            metallic = True
        if a_type == 'B':
            boron_count += 1
        if a_type in ['Cu', 'Fe', 'Mn', 'Ni']:
            ok_count += 1
        if a_type == 'C':
            c_count += 1

    # metallic or non metallic with more than 7 boron atoms will be excluded
    # metallic molecules with 0 C atoms and cu,fe,mn, or ni won't be excluded
    # not RDKit sanitizable will be excluded
    if (metallic and not (c_count == 0 and ok_count != 0)) or (not metallic and boron_count > 7) or rdkit_fails:
        exclude = True

    return exclude
