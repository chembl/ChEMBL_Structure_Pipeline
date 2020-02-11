""" The ChEMBL structurepipeline

# Checking structures for (potential) problems.

The primary checking function here is check_molblock(), which takes a Mol block as 
an argument and returns a tuple of (penalty, description tuples) describing what
was found:

    >>> check_molblock(''' 
    ...   Mrv1810 02151908462D           
    ...  
    ...   4  3  0  0  0  0            999 V2000 
    ...     2.2321    4.4196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 
    ...     3.0023    4.7153    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 
    ...     1.4117    4.5059    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0 
    ...     1.9568    3.6420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 
    ...   1  2  1  1  0  0  0 
    ...   1  3  1  0  0  0  0 
    ...   1  4  1  0  0  0  0 
    ... M  END 
    ... ''')                                                                                                            
    ((5, 'InChi_RDKit/Mol stereo mismatch'),)
    >>>

"Clean" Mol blocks return an empty tuple:

    >>> check_molblock(''' 
    ...   Mrv1810 02151908462D           
    ...  
    ...   4  3  0  0  0  0            999 V2000 
    ...     2.2321    4.4196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 
    ...     3.0023    4.7153    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0 
    ...     1.4117    4.5059    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0 
    ...     1.9568    3.6420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 
    ...   1  2  1  1  0  0  0 
    ...   1  3  1  0  0  0  0 
    ...   1  4  1  0  0  0  0 
    ... M  END 
    ... ''')                                                                                                            
    ()
    >>>

# Standardizing structures:

    >>> mb = standardize_molblock('''
    ...   Mrv1810 07121910172D          
    ... 
    ...   4  3  0  0  0  0            999 V2000
    ...    -2.5038    0.4060    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0
    ...    -2.5038    1.2310    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
    ...    -3.2182   -0.0065    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    ...    -1.7893   -0.0065    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    ...   1  2  1  0  0  0  0
    ...   1  3  1  0  0  0  0
    ...   1  4  1  4  0  0  0
    ... M  CHG  2   2  -1   3   1
    ... M  END
    ... ''')
    >>> print(mb)
    <BLANKLINE>
         RDKit          2D
    <BLANKLINE>
      4  3  0  0  0  0  0  0  0  0999 V2000
       -2.5038    0.4060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -2.5038    1.2310    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
       -3.2182   -0.0065    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
       -1.7893   -0.0065    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0
      1  3  1  0
      1  4  1  0
    M  END
    <BLANKLINE>
   >>>

# Standardization + salt stripping :
    >>> mb,exclude = get_parent_molblock('''
    ...   Mrv1810 07121910262D          
    ... 
    ...   3  1  0  0  0  0            999 V2000
    ...    -5.2331    1.1053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    ...    -4.5186    1.5178    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    ...    -2.8647    1.5789    0.0000 Cl  0  5  0  0  0  0  0  0  0  0  0  0
    ...   1  2  1  0  0  0  0
    ... M  CHG  2   2   1   3  -1
    ... M  END
    ... ''')
    >>> print(mb)
    <BLANKLINE>
         RDKit          2D
    <BLANKLINE>
      2  1  0  0  0  0  0  0  0  0999 V2000
       -5.2331    1.1053    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -4.5186    1.5178    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0
    M  END
    <BLANKLINE>
    >>>

"""
from .checker import check_molblock
from .standardizer import standardize_molblock, standardize_mol
from .standardizer import get_parent_molblock, get_parent_mol
__version__ = "1.0.0"

#
#  Copyright (c) 2019 Greg Landrum
#  All rights reserved.
#
#  This file is part of the ChEMBL_StructurePipeline project.
#  The contents are covered by the terms of the MIT license
#  which is included in the file LICENSE, found at the root
#  of the source tree.
