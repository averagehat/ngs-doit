from typing_extensions import Protocol, TypedDict
from typing import TypeVar, Generic, List, Tuple, Union
# doesn't really work becuase the semantics (what we really mean)
# is that targets is a list which includes at least one Bam file. 
import pathlib


T = TypeVar('T', bound='FPath')
class FPath(pathlib.PurePosixPath):

    def add_suffix(self: 'T', suf: str) -> 'T':
        return self.with_suffix(self.suffix + suf)

    @property
    def fasta(self) -> 'Fasta':
        return Fasta(str(self))

    @property
    def fastq(self) -> 'Fastq': 
        return Fastq(self.add_suffix('.fastq'))

    @property
    def bam(self) -> 'Bam': 
        return Bam(self.add_suffix('.bam'))

    @property
    def bai(self) -> 'Bai': 
        return Bai(self.add_suffix('.bai'))

    @property
    def bai(self) -> 'Sam': 
        return Sam(self.add_suffix('.sam'))

    @property
    def vcf(self) -> 'VCF': 
        return VCF(self.add_suffix('.vcf'))

    @property
    def png(self) -> 'PNG': 
        return PNG(self.add_suffix('.png'))
        # will replace old suffix
        # requires that you start suffix with '.'
        return Fastq(self.with_suffix('.fastq'))

class Fastq(FPath): ...
class Bam(FPath): ...
class Sam(FPath): ...
class Bai(FPath): ...
class VCF(FPath): ...
class Fasta(FPath): ...
class PNG(FPath): ...

PEJob = TypedDict('PEJob', { 'dependencies' : List[Tuple[Fastq, Fastq]] })
PEMap = TypedDict('PEMap', { 'targets' : List[Bam], 'dependencies' : List[Tuple[Fastq, Fastq]] }) 
