from rdkit import RDLogger
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit import Chem
import argparse
import warnings
import sys
warnings.simplefilter('ignore', category=RuntimeWarning)


class RunningException(Exception):
    pass


class MolFromSmiles(RunningException):
    pass


class HasPains(RunningException):
    pass


class CouldNotBeStandartized(RunningException):
    pass


class Filter:
    p = FilterCatalogParams.FilterCatalogs.PAINS
    A = FilterCatalogParams.FilterCatalogs.PAINS_A
    B = FilterCatalogParams.FilterCatalogs.PAINS_B
    C = FilterCatalogParams.FilterCatalogs.PAINS_C


def main():
    parser = argparse.ArgumentParser(
        description='Sanitize smiles using chembl_structure_pipeline and RDKit PAINS filters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-s', '--standartize', action='store_true', default=True,
                        help='Whether to perform standartization of input SMILES')
    parser.add_argument('-p', action='store_true', default=True,
                        help='Filter molecules using all PAINS filters together')
    parser.add_argument('-A', action='store_true', default=False,
                        help='Filter molecules using all PAINS_A filter separately')
    parser.add_argument('-B', action='store_true', default=False,
                        help='Filter molecules using all PBINS_B filter separately')
    parser.add_argument('-C', action='store_true', default=False,
                        help='Filter molecules using all PCINS_C filter separately')
    parser.add_argument('input', metavar='INPUT',
                        help='Input file (with SMILES as first column)')
    parser.add_argument('--strict', action='store_true', default=False,
                        help='Whether to raise an exception on first error')
    parser.add_argument('--header', action='store_true', default=False,
                        help='Indicate that the input file contains header')
    parser.add_argument('--verbose', action='store_true', default=False,
                        help='Whether to print all RDKit warnings to stdout')
    parser.add_argument('--stderr', action='store_true', default=False,
                        help='Whether to print filtered molecules to stderr')
    args = parser.parse_args()
    if not args.verbose:
        from rdkit.rdBase import BlockLogs
        block = BlockLogs()
        RDLogger.DisableLog('rdApp.*')

    params = FilterCatalogParams()
    if args.p:
        params.AddCatalog(Filter.p)
    if args.A:
        params.AddCatalog(Filter.A)
    if args.B:
        params.AddCatalog(Filter.B)
    if args.C:
        params.AddCatalog(Filter.C)
    catalog = FilterCatalog(params)

    def has_pains(mol):
        entry = catalog.GetFirstMatch(mol)
        return bool(entry), entry  # will be empty if no match found

    with open(args.input, 'r') as fin:
        from chembl_structure_pipeline import standardize_mol

        if args.header:
            line = next(fin)
            print(line, end='')

        for idx, line in enumerate(fin):
            smiles = line.strip().split(maxsplit=1)[0]
            try:
                mol = Chem.MolFromSmiles(smiles)
                if not mol:
                    raise MolFromSmiles(
                        'Line {}, SMILES: {}'.format(idx, smiles))

                mol_std = standardize_mol(mol)
                if not mol_std:
                    raise CouldNotBeStandartized(
                        'Line {}, SMILES: {}'.format(idx, smiles))

                result, entry = has_pains(mol_std)
                if result:
                    group = entry.GetProp('Scope')
                    raise HasPains(
                        'Line {}, SMILES: {}, type: {}'.format(idx, smiles, group))

                print(line, end='')

            except RunningException:
                if args.strict:
                    raise
                if args.stderr:
                    print(line, end='', file=sys.stderr)


if __name__ == '__main__':
    main()
