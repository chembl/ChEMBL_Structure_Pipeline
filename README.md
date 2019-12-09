[![Build Status](https://travis-ci.org/chembl/ChEMBL_Structure_Pipeline.svg?branch=master)](https://travis-ci.org/chembl/ChEMBL_Structure_Pipeline)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# ChEMBL Structure Pipeline

ChEMBL protocols used to standardise and salt strip molecules. First used in ChEMBL 26.

Check the [wiki](https://github.com/chembl/ChEMBL_Structure_Pipeline/wiki) for some implementation details.

## Installation (it requires RDKit to work)

From source:

    git clone https://github.com/chembl/ChEMBL_Structure_Pipeline.git
    pip install ./ChEMBL_Structure_Pipeline

From conda:

```bash
conda install -c chembl structurepipeline
```

## Usage

### Standardise molblock


```python
from structurepipeline import standardizer

o_molblock = """
  SciTegic07111214002D
 22 19  0  0  0  0            999 V2000
   12.9184   -8.1544    0.0000 C   0  0
   13.6343   -7.7410    0.0000 N   0  0
   12.2024   -7.7410    0.0000 O   0  0
   12.9184   -8.9811    0.0000 N   0  0
   14.3503   -8.1544    0.0000 C   0  0
   15.0663   -7.7410    0.0000 C   0  0
   15.7784   -8.1460    0.0000 C   0  0
   16.4924   -7.7327    0.0000 C   0  0
   17.2064   -8.1418    0.0000 C   0  0
   17.9204   -7.7285    0.0000 C   0  0
   18.6343   -8.1376    0.0000 C   0  0
   19.3483   -7.7243    0.0000 C   0  0
   20.0623   -8.1335    0.0000 C   0  0
   20.7763   -7.7202    0.0000 C   0  0
   21.4902   -8.1293    0.0000 C   0  0
   22.2042   -7.7160    0.0000 C   0  0
   22.9207   -8.1284    0.0000 N   0  0
   23.6362   -7.7143    0.0000 C   0  0
   24.3526   -8.1267    0.0000 O   0  0
   23.6352   -6.8875    0.0000 N   0  0
   14.3464   -5.9039    0.0000 Cl  0  0
   24.1332   -6.2672    0.0000 Cl  0  0
  2  5  1  0
 11 12  1  0
  1  2  1  0
 12 13  1  0
  5  6  1  0
 13 14  1  0
 14 15  1  0
  6  7  1  0
 15 16  1  0
  1  3  1  0
 16 17  1  0
  7  8  1  0
 17 18  1  0
 18 19  1  0
  8  9  1  0
 18 20  2  0
  1  4  2  0
  9 10  1  0
 10 11  1  0
M  END"""

std_molblock = standardizer.standardize_molblock(o_molblock)
```

### Get the parent molecule


```python
from structurepipeline import standardizer


o_molblock = """
  SciTegic06221110232D

 23 22  0  0  0  0            999 V2000
    3.5724    1.2375    0.0000 C   0  0
    2.8579    0.8250    0.0000 C   0  0
    2.8579    0.0000    0.0000 O   0  0
    2.1434    1.2375    0.0000 O   0  0
    1.4289    0.8250    0.0000 C   0  0
    1.4289    0.0000    0.0000 C   0  0
    0.7145   -0.4125    0.0000 C   0  0
    0.0000    0.0000    0.0000 C   0  0
    0.0000    0.8250    0.0000 C   0  0
    0.7145    1.2375    0.0000 C   0  0
    0.7145    2.0625    0.0000 C   0  0
    1.4289    2.4750    0.0000 O   0  0
    0.0000    2.4750    0.0000 O   0  0
    2.4973   -3.1821    0.0000 O   0  0
    1.7828   -2.7696    0.0000 C   0  0
    1.0683   -3.1821    0.0000 C   0  0
    0.3539   -2.7696    0.0000 C   0  0
   -0.3606   -3.1821    0.0000 C   0  0
   -1.0751   -2.7696    0.0000 C   0  0
   -1.7895   -3.1821    0.0000 C   0  0
   -2.5040   -2.7696    0.0000 N   0  0
    1.7828   -1.9446    0.0000 O   0  0
    1.0683   -4.0071    0.0000 N   0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  2  0
  5 10  1  0
 10 11  1  0
 11 12  1  0
 11 13  2  0
 14 15  1  0
 15 16  1  0
 16 17  1  0
 17 18  1  0
 18 19  1  0
 19 20  1  0
 20 21  1  0
 15 22  2  0
 16 23  1  0
M  END
"""

parent_molblock, _ = standardizer.get_parent_molblock(o_molblock)
```
