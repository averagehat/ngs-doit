from Bio.SeqIO.QualityIO import FastqGeneralIterator as FQGen
from itertools import chain
from functools import partial
from typing import NamedTuple, Iterator, Optional, Tuple
from typing_extensions import Final
class FqPair(NamedTuple):
    fwd: str
    rev: str

# simplefastq parser (FastqGeneralIterator) returns this type
# form is (header, sequence, quals)
Rec = Tuple[str, str, str]
PHRED_OFFSET: Final = 33

def fastq_filter(rps: FqPair, indexes: Optional[FqPair], outs: FqPair, minIndexBQ: int, keepNs: bool) -> None: 

   ORD_MIN = minIndexBQ + PHRED_OFFSET

   def keep_rec(r1: Rec, r2: Rec, i1: Optional[Rec], i2: Optional[Rec]) -> Optional[Tuple[Rec, Rec]]:
       if i1 and i2: 
           quals = chain(i1[2], i2[2])
           uni_codes = map(ord, quals)
           qual_good = all( x >= ORD_MIN for x in uni_codes )
       bases = chain(r1[1], r2[1])
       keep = qual_good and (keepNs or (not ('N' in bases)))
       return (r1, r2) if keep else None

   with open(outs.fwd, 'w') as fwdout, \
           open(outs.rev, 'w') as revout, \
           open(rps.fwd) as fwdin, \
           open(rps.rev) as revin:

       if indexes:
           with open(indexes.fwd) as idxfwd, open(indexes.rev) as idxrev:
                result = map(keep_rec, FQGen(fwdin), FQGen(revin), FQGen(idxfwd), FQGen(idxrev))
       else:
            func = partial(keep_rec, i1=None, i2=None)
            result = map(func, FQGen(fwdin), FQGen(revin))

       kept = filter(None, result)
       for fwd, rev in kept:
           fwdout.write("@%s\n%s\n+\n%s\n" % fwd )
           revout.write("@%s\n%s\n+\n%s\n" % rev )
