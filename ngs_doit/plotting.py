# type: ignore
#_, poss, depths = itertools.tee(csv.DictReader(open('depth.hist'),  fieldnames=('chrom', 'pos', 'depth')))
#print(list(csv.reader(open('depth.hist'), delimiter='\t'))[:4])
##_, poss, depths = itertools.tee(csv.reader(open('depth.hist'), delimiter='\t'))
#reader = csv.reader(open('depth.hist'), delimiter='\t')
#_, depthX, depthY = zip(*list(reader))
#rows = csv.DictReader(open('qual.hist'), delimiter='\t')
##qualX, qualY = (x[0], avg(x[1] + x[3])
#qualX, qualY = zip(*[(x['#BaseNum'], avg(x['Read1_linear'], x['Read2_linear'] ) ) for x in rows])
#_,  xs, ys = itertools.tee(reader, 3)
#xs, ys = (int(x[1]) for x in xs), (int(y[2]) for y in ys)
import pysam 
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import click 

import os

# just types
from plotly.graph_objs import Figure 
from pysam.libcalignedsegment import PileupColumn
from pysam.libcalignmentfile import AlignmentFile
from dataclasses import dataclass
from typing import Sequence, Iterator, Union, Dict
from typing_extensions import Literal

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
    #fig.update_layout(yaxis_range=(0, 100))
    return fig

def save_fig(fig: Figure, basepath: str, ft: FigType) -> None:
    fn = f"{basepath}.{ft}" #NOTE: this will replace the last suffix.
    with open(fn, 'w') as out:
        if ft is 'html':
            fig.write_html(out)
        elif ft is 'png':
            fig.write_image(out)
        elif ft is 'json':
            fig.write_json(out)
        else:
            assert "this should never happen."

MAX_BASE_QUALITY = 50

@click.command()
@click.argument('bampath', type=click.Path(exists=True, readable=True))
@click.argument('outpath', type=click.Path(exists=False, writable=True, dir_okay=True))
@click.option('--format', multiple=True, type=click.Choice(['html', 'png', 'json']), default=['html', 'json'])
@click.option('--orphans/--no-orphans', is_flag=True, default=True)
@click.option('--overlaps/--no-overlaps', is_flag=True, default=False)
@click.option('--minbq', type=click.IntRange(0, MAX_BASE_QUALITY), default=0)
def plot_coverage(bampath: str, outpath: str, format: Sequence[FigType], orphans: bool, overlaps: bool, minbq: int) -> None:
    assert os.path.exists(bampath +'.bai'), f"Expedcted index file {bampath + '.bai'} not found."
    opts = PileupOptions(minBaseQual=minbq, ignoreOrphans = not orphans, ignoreOverlaps = not overlaps)
    cvgs = full_coverage(bampath, opts)
    fig = coverage_plot(cvgs) 
    for ft in format:
        save_fig(fig, outpath, ft)

if __name__ == '__main__':
    plot_coverage()

