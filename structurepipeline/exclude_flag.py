from rdkit import Chem


# Zn not in the list as we have some compounds in ChEMBL
# most of them are simple salts
METAL_LIST = [
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Ga", "Y", "Zr", "Nb", "Mo",
    "Tc", "Ru", "Rh", "Pd", "Cd", "In", "Sn", "La", "Hf", "Ta", "W", "Re", "Os", "Ir",
    "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Ac", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", 
    "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Th", "Pa", "U", "Np", "Pu", "Am",
    "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr","Ge", "Sb",
]


def exclude_flag(molfile, includeRDKitSanitization=True):
    """
    Rules to exclude structures. Ported and updated from original PP protocol.
    """
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

    for atom in mol.GetAtoms():
        a_type = atom.GetSymbol()
        if a_type in METAL_LIST:
            metallic = True
        if a_type == "B":
            boron_count += 1
        if a_type in ["Cu", "Fe", "Mn", "Ni"]:
            ok_count += 1
        if a_type == "C":
            c_count += 1

    # metallic or non metallic with more than 7 boron atoms will be excluded
    #  - problem of depicting borane compounds
    # RDKit non sanitizable molecules will also be excluded
    # metallic molecules with no C atoms and Cu, Fe, Mn, or Ni won't be excluded
    if (
        (metallic and not (c_count == 0 and ok_count != 0))
        or (not metallic and boron_count > 7)
        or rdkit_fails
    ):
        exclude = True

    return exclude
