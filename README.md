[![CI Testing](https://github.com/chembl/ChEMBL_Structure_Pipeline/workflows/CI/badge.svg)](https://github.com/chembl/ChEMBL_Structure_Pipeline/actions?query=workflow%3ACI+branch%3Amain)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# ChEMBL Structure Pipeline

ChEMBL protocols used to standardise and salt strip molecules. First used in ChEMBL 26.

Check the [wiki](https://github.com/chembl/ChEMBL_Structure_Pipeline/wiki) and paper[[1]](#1) for a detailed description of the different processes.

## Installation (it requires RDKit to work)

From source:

    git clone https://github.com/chembl/ChEMBL_Structure_Pipeline.git
    pip install ./ChEMBL_Structure_Pipeline

Using conda:

```bash
conda install -c conda-forge chembl_structure_pipeline
```

## Usage

### Standardise a compound [(info)](https://github.com/chembl/ChEMBL_Structure_Pipeline/wiki/Work-done-by-each-step#standardize_molblock)


```python
from chembl_structure_pipeline import standardizer

o_molblock = """
  Mrv1810 07121910172D          

  4  3  0  0  0  0            999 V2000
   -2.5038    0.4060    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
   -2.5038    1.2310    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
   -3.2182   -0.0065    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   -1.7893   -0.0065    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  4  0  0  0
M  CHG  2   2  -1   3   1
M  END
"""

std_molblock = standardizer.standardize_molblock(o_molblock)
```

### Get the parent compound [(info)](https://github.com/chembl/ChEMBL_Structure_Pipeline/wiki/Work-done-by-each-step#get_parent_molblock)


```python
from chembl_structure_pipeline import standardizer

o_molblock = """
  Mrv1810 07121910262D          

  3  1  0  0  0  0            999 V2000
   -5.2331    1.1053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.5186    1.5178    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   -2.8647    1.5789    0.0000 Cl  0  5  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  CHG  2   2   1   3  -1
M  END
"""

parent_molblock, _ = standardizer.get_parent_molblock(o_molblock)
```

### Check a compound [(info)](https://github.com/chembl/ChEMBL_Structure_Pipeline/wiki/Work-done-by-each-step#checkmolecule)

The checker assesses the quality of a structure. It highlights specific features or issues in the structure that may need to be revised. Together with the description of the issue, the checker process returns a penalty score (between 0-9) which reflects the seriousness of the issue (the higher the score, the more critical is the issue)

```python
from chembl_structure_pipeline import checker

o_molblock = """ 
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

issues = checker.check_molblock(o_molblock)
```

## References
<a id="1">[1]</a> 
Bento, A.P., Hersey, A., Félix, E. et al. An open source chemical structure curation pipeline using RDKit. J Cheminform 12, 51 (2020). https://doi.org/10.1186/s13321-020-00456-1

