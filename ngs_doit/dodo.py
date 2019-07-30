from types import Fasta, Fastq, Bam, Bai, VCF, PNG ,FPath, Sam
from typing import List, Tuple, Dict, Callable,Union, Any,Optional
from typing_extensions import TypedDict


Job = Dict[str, Any]

# the type for get_var is something like
# def get_var(s: str, default: T) -> Union[str, T] (the same as dict.get())
from doit import get_var # type: ignore     

import pathlib
import subprocess
import os

import yaml

# could use "Old" and f"" style string will work 
# TODO: the below doesn't work@
DOIT_CONFIG = {'action_string_formatting': 'both'}

Args = TypedDict('Args', 
        { 'config' : str,
          'outdir' : str,
          'ref'    : Optional[str],
          'R1'     : Optional[str],
          'R2'     : Optional[str] })
args: Args = {"config": get_var('config', 'config.yaml'),
          "outdir": get_var('outdir', './re-out'),
          "ref"   : get_var('ref', None),
          "R1"    : get_var('r1', None),
          'R2'    : get_var('r2', None)}


# these file names could be dynamically named after the sample by adding a "sample" value to
# the `Args` dict and using that within each relevant task's closure! 
[trimmed_r1, trimmed_r2, up_trimmed1, up_trimmed2] = \
        [ FPath('trimmed').fastq, FPath('Trimmed.R2').fastq, \
        Fastq('unpaired_fw.fq'), Fastq('unpaired_rev.fq')]
ref = Fasta(args['ref']) # type: ignore
index_files: List[FPath] = [ref.add_suffix(suf) for suf in  { '.amb', '.sa', '.ann', '.bwt', '.pac'} ]
unsorted_bam = Bam('bwa.bam')
sorted_bam = Bam("sorted.bam")
index_bai = sorted_bam.bai
# could name below in a standard way but this is better

def task_trimmomatic() -> Union[Dict[str, Any], bool]:
    # could also have this at the "main level" of the file :(
    cfg = yaml.load(open(args['config']))
    q  = cfg['trim_reads']['q']['default']
    hc = cfg['trim_reads']['headcrop']['default']
    r1 = args['R1']
    r2 = args['R2'] 
    if r1 is None or r2 is None:
        return False # some other way of indicating failure

    def run_trim(targets: List[str], dependencies: List[str]) -> None:
        subprocess.getoutput(f"trimmomatic PE {' '.join(dependencies)} \
                {targets[0]} unpaired_fw.fq {targets[1]} unpaired_rev.fq \
                LEADING:{q} TRAILING:{q} HEADCROP:{hc}")

    return { 'targets' : [ FPath('trimmed').fastq, FPath('Trimmed.R2').fastq],
             'file_dep' : [r2, r1],
             'actions' : [run_trim] }

def task_bwa_index() -> Job:

    return { 'targets' : index_files,
             'file_dep' : [ref] ,
             'actions' : [ "bwa index %(dependencies)s" ] 
             }

def task_bwa_mem() -> Job:
    # needs to dynamically whether or not unpaired_fwd from trimmed reads are empty files or not; alternatively, whether or not they exist.
    # could just have merged-bam be the ending output, and so no worries. problem with this is that paired and unpaired could be paralleizable but meh

    # these are actually paths right?
    # let sorting be a seperator process 
    def bwa_mem(targets: List[str], dependencies: List[str]) -> None:
        nonempty: Callable[[str],bool] = lambda f: os.stat(f).st_size > 0 
        r1, r2, u1, u2 = targets
        paired_bam =   Bam('paired.bam')
        unpaired_bam = Bam('unpaired.bam')
        final_bam = Bam(dependencies[0]) # eh? == unsorted_bam
    # could use overload types somehow
        def map_paired(outpath: Bam) -> None:
            subprocess.getoutput(f"bwa mem {r1} {r2} -o {outpath}")
        if nonempty(dependencies[0]) or nonempty(dependencies[1]):
            unpaired = Fastq("unpaired.fastq")
            subprocess.getoutput(f"cat {u1} {u2} > {unpaired}")
            map_paired(paired_bam)
            subprocess.getoutput(f"samtools merge {unpaired_bam} {paired_bam} > {final_bam}")
            os.remove(unpaired_bam)
            os.remove(paired_bam) # otherwise the targets change!
        else:
            map_paired(final_bam) # but doesn't sort :(

    return { 'targets' : unsorted_bam,
             'file_dep' : index_files + [trimmed_r1, trimmed_r2, up_trimmed1, up_trimmed2], 
             'actions' : [ bwa_mem ] }

def task_sort_bam() -> Job:
    # how can these targets and dependencies be parameterizable like they are in, say, `dagr`?
    return { 'targets' : sorted_bam,
            'file_dep' :  unsorted_bam,
            'actions'  : [ "samtools sort %(dependencies)s -f %(targets)s" ] }

def task_index_bam() -> Job: 

    return { 'targets'  : index_bai,
             'file_dep' : sorted_bam,
             'actions'  : [ "samtools index %(dependencies)s" ] }


def task_fqstats() -> Job:
    cfg = yaml.load(open(args['config']))
    skip_stats = cfg['fqstats']['skip']['default'] # or accept commandline arg

    if skip_stats:
        return { 'targets' : [], 'file_dep' : [], 'actions': []}
    else: 
        return { 'targets' : ["foo.fqstat"], 'file_dep' : [], 'actions' : ['fqstats bleh bleh'] }


