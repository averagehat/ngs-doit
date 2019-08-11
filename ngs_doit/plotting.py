# type: ignore
import pysam 
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import click 
from Bio.SeqIO.QualityIO import FastqGeneralIterator

import os
from operator import itemgetter as atindex
from collections import Counter
from itertools import chain
from pathlib import PurePosixPath

# just types
from ngs_doit.tools import Rec 
from typing import List, IO
from plotly.graph_objects import Figure 
from pysam.libcalignedsegment import PileupColumn
from pysam.libcalignmentfile import AlignmentFile
from dataclasses import dataclass
from typing import Sequence, Iterator, Union
from typing_extensions import Literal, Final
from ngs_doit.custom_types import FileDeps, Targets

FigType = Literal['html', 'png', 'json']
Ref = str


@dataclass
class Coverage:
    depth: int
    avgBaseQual: float
    referencePos: int
    reference: Ref
    # ref field?

@dataclass
class PileupOptions:
    # minDepth: int
    minBaseQual: int
    ignoreOrphans: bool
    ignoreOverlaps: bool

def pileup(bam: AlignmentFile, ref: Ref, start: int, stop: int, opts: PileupOptions) -> Iterator[PileupColumn]:
    assert ref in bam.references
    assert start >= 0 and stop <= bam.get_reference_length(ref) 
    assert bam.check_index()
    assert not bam.closed
    assert bam.is_bam

    BIG_INT = 999_999_999

    for x in bam.pileup(ref, start, stop, max_depth=BIG_INT,
               min_base_quality=opts.minBaseQual,
               ignore_orphans=opts.ignoreOrphans,
               ignore_overlaps=opts.ignoreOverlaps):
        yield x 

def avg(xs: Sequence[Union[int, float]]) -> float:
    return sum(xs) / len(xs)

def col2coverage(x: PileupColumn) -> Coverage:
    return Coverage(referencePos=x.reference_pos,
                    depth=x.nsegments,
                    avgBaseQual=avg(x.get_query_qualities()),
                    reference=x.reference_name)

def full_coverage(bampath: str, opts: PileupOptions) -> Sequence[Coverage]:
    with pysam.AlignmentFile(bampath) as bam:
        for ref in bam.references:
            pileupgen = pileup(bam, ref, 0, bam.get_reference_length(ref), opts)
            cvgs = map(col2coverage, pileupgen)
            for cvg in cvgs:
                yield cvg

def coverage_plot(cs: Sequence[Coverage]) -> Figure: 
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    xs, depths, quals = zip(*[(x.referencePos, x.depth, x.avgBaseQual) for x in cs])
    fig.add_trace(go.Scatter(x=xs, y=depths, fill='tozeroy', stackgroup='one', 
        mode='lines',
        line=dict(width=0.5, color='rgb(131, 90, 241)'))) # fill down to xaxis
    fig.add_trace(go.Scatter(x=xs, y=quals, fill='tonexty', stackgroup='one',
        mode='lines',
        line=dict(width=0.5, color='rgb(111, 231, 219)')), secondary_y=True) # fill down to xaxis
    return fig

def save_fig(fig: Figure, basepath: str, ft: FigType) -> None:
    #fn = f"{basepath}.{ft}"
    out = f"{basepath}.{ft}"
    if ft == 'html':
        fig.write_html(out)
    elif ft == 'png':
        fig.write_image(out)
    elif ft == 'json':
        fig.write_json(out)
    else:
        assert False, "this should never happen."



def fastq_lengths(f: IO[str]) -> Iterator[int]:
       recs: Rec = FastqGeneralIterator(f)
       lengths = map(len, map(atindex(1), recs))
       return lengths

def plot_length_dist(targets: List[str], dependencies: List[str]) -> None:
# could calculate any distribution here, i.e. avg quality
   #with open(fqpath) as f:
   lengths = chain.from_iterable(map(fastq_lengths, map(open, dependencies)))
   counts = Counter(lengths)
   xs = sorted(counts.keys())
   ys = [counts[k] for k in xs]
   fig = go.Figure(go.Scatter(x=xs,y=ys))
   for p in map(PurePosixPath, targets):
       save_fig(fig, p.with_suffix(''), p.suffix[1:])
           
MAX_BASE_QUALITY: Final = 41
#@click.command() # type: ignore
#@click.argument('bampath', type=click.Path(exists=True, readable=True)) # type: ignore
#@click.option('--out', type=click.Path(exists=False, writable=True, dir_okay=True), multiple=True) # type: ignore
#@click.option('--orphans/--no-orphans', is_flag=True, default=True) # type: ignore
#@click.option('--overlaps/--no-overlaps', is_flag=True, default=False) # type: ignore
#@click.option('--minbq', type=click.IntRange(0, MAX_BASE_QUALITY), default=0) # type: ignore
#def plot_coverage(targets: Targets, dependencies: FileDeps, orphans: bool, overlaps: bool, minbq: int) -> None:
def plot_coverage(bampath: str, targets: List[str], orphans: bool, overlaps: bool, minbq: int) -> None:
    #bampath = dependencies[0]
    assert os.path.exists(bampath +'.bai'), f"Expedcted index file {bampath + '.bai'} not found."
    opts = PileupOptions(minBaseQual=minbq, ignoreOrphans = not orphans, ignoreOverlaps = not overlaps)
    cvgs = full_coverage(bampath, opts)
    fig = coverage_plot(cvgs) 
    for p in map(PurePosixPath, targets):
        save_fig(fig, p.with_suffix(''), p.suffix[1:])

if __name__ == '__main__':
    plot_coverage()

