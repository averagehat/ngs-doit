import re
import os
from pathlib import Path, PurePosixPath
from typing import Callable, Tuple, Sequence, Optional, Dict, Any, Union, List, NoReturn, Iterator, IO
from typing_extensions import TypedDict
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from itertools import chain
#from ngs_doint import Job

def get_index(fn: Path) -> Optional[Path]:
    ''' returns index path (ie _I1_ or _I2_)  or none.'''
    index = re.sub(r'_R([12])_', r'_I\1_', str(fn))
    return Path(index) if os.path.exists(index) else None

def task_ngs_filter() -> Any: #Job: 
    filtered1_fq, filtered2_fq = 'filtered.r1.fq', 'filtered.r2.fq'
    input1_fq, input2_fq = 'r1.f1', 'f2.fq' # args['R1'], args['R2'] 
    def filter_action(targets: List[str], dependencies: List[str]) -> bool:
        # index quality min (avg? I don't remember)
        # avg qual min
        # N maximum count
        return False
    return { 'targets' : [filtered1_fq, filtered2_fq],
             'file_dep' : [input1_fq, input2_fq],
             'actions' : [filter_action] }

# below is the type returned by fastqgeneraliterator
Rec = Tuple[str, str, str]
CHUNKSIZE = 10000   # Set this to whatever you feel reasonable

def fq_filter_parallel(r1p, i1p, r2p, i2p, outR1, outR2, minIndexBQ, keepNs):
   togen = lambda x: FastqGeneralIterator(open(x))
   records = zip(togen(r1p), togen(i1p), togen(r2p), togen(i2p))
   with open(outR1, 'w') as r1File, \
           open(outR2, 'w') as r2File:

       ORD_MIN = minIndexBQ + 33

       def keep_rec(x: Tuple[Rec, Rec, Rec, Rec]) -> bool:
           Iqual1, Iqual2 = x[1][2], x[3][2]
           bases1, bases2 = x[0][1], x[2][1]
           qual_good = map(lambda x: x >= ORD_MIN, \
                       map(ord, chain(Iqual1, Iqual2)))
           result = (keepNs or (not ('N' in chain(bases1, bases2)))) \
                   and all(qual_good)
           return result

       good = filter(keep_rec, records)

       for r1, _, r2, _ in good:
           r1File.write("@%s\n%s\n+\n%s\n" % r1 )
           r2File.write("@%s\n%s\n+\n%s\n" % r2 )
