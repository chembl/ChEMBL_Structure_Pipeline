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

"""
from .checker import check_molblock
__version__ = "0.0.1"

#
#  Copyright (c) 2019 Greg Landrum
#  All rights reserved.
#
#  This file is part of the ChEMBL_StructurePipeline project.
#  The contents are covered by the terms of the MIT license
#  which is included in the file LICENSE, found at the root
#  of the source tree.
