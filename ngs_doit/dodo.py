from ngs_doit.pipetypes import Fasta, Fastq, Bam, Bai, VCF, PNG ,FPath, Sam, P
from typing import List, Tuple, Dict, Callable,Union, Any,Optional, TypeVar, Type #, Set
from typing_extensions import TypedDict, Literal

# the type for get_var is something like
# def get_var(s: str, default: T) -> Union[str, T] (the same as dict.get())
from doit import get_var # type: ignore     

import pathlib
import subprocess
import os

import yaml

#Job = Dict[str, Any]
#P = TypeVar('P', bound=FPath)
#Action = Union[str, Callable[[List[Type[P]], List[FPath]], None]]
Action = Union[str, Callable[[List[FPath], List[FPath]], ActionReturn]]
ActionReturn = Union[Failure, None, Dict[str, bool]]

# could use "Old" and f"" style string will work 
# TODO: the below doesn't work@
DOIT_CONFIG = {'action_string_formatting': 'both'}

Args = TypedDict('Args', 
        { 'config' : str,
          'outdir' : str,
          'ref'    : Optional[str],
          'R1'     : Optional[str],
          'R2'     : Optional[str],
          'sample' : Optional[str]})
args: Args = {"config": get_var('config', 'config.yaml'),
          "outdir": get_var('outdir', './re-out'),
          "ref"   : get_var('ref', None),
          "R1"    : get_var('r1', None),
          'R2'    : get_var('r2', None),
          'sample': get_var('sample', None)}

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
IntValue = TypedDict('IntValue', { 'default' : int })
StringValue = TypedDict('StringValue', { 'default' : str })
BoolValue = TypedDict('BoolValue', { 'default' : bool })
Trim = TypedDict('Trim', {'q' : IntValue , 'headcrop' : IntValue })
Fqstats = TypedDict('Fqstats', {'skip' : BoolValue })
Config = TypedDict('Config', 
        { 'trim_reads' :  Trim,
    'fqstats' : Fqstats
    })
from typing import cast

Job = TypedDict('Job', 
        { 'file_dep' : List[FPath],
          'targets' : List[FPath],
          'actions' : List[Action] })

GetArgs = Dict[str, Tuple[str, str]]
# could do some cute b.s. with kwarg typing & generics but generics unsupported for typeddict
BoolArgJob = TypedDict('BoolArgJob', { 'getargs' : GetArgs, 
         'file_dep' : List[FPath],
          'targets' : List[FPath],
          'actions' : List[GAction[bool]] })

class CmdAction: ...
from typing import TypeVar
T = TypeVar('T')
# from doit.action import CmdAction  # type:ignore
GAction = Union[CmdAction, Callable[[List[FPath], List[FPath], T], ActionReturn]] 
Failure = Literal[False]

def task_trimmomatic() -> Union[Job, Failure]:
    # could also have this at the "main level" of the file :(
    cfg: Config = yaml.load(open(args['config']))
    q  = cfg['trim_reads']['q']['default']
    hc = cfg['trim_reads']['headcrop']['default']
    r1_ = args['R1']
    r2_ = args['R2'] 
    if r1_ is None or r2_ is None:
        return False # some other way of indicating failure
    r1, r2 = Fastq(r1_), Fastq(r2_)
    def run_trim(targets: List[FPath], dependencies: List[FPath]) -> Dict[str, bool]:
        nonempty: Callable[[FPath],bool] = lambda f: os.stat(f).st_size > 0 
        subprocess.getoutput(f"trimmomatic PE {dependencies[0]} {dependencies[1]} \
                %(targets)s LEADING:{q} TRAILING:{q} HEADCROP:{hc}")
        return { 'unpaired_reads' : nonempty(targets[1]) or nonempty(targets[3]) }

    return { 'targets' : [ FPath('trimmed.R1').fastq,  Fastq('unpaired_trimmed_R1.fastq'), FPath('Trimmed.R2').fastq, Fastq('unpaired_trimmed_R2.fastq')],
             'file_dep' : [r2, r1],
             'actions' : [run_trim] }


def task_bwa_index() -> Job:

    return { 'targets' : index_files,
             'file_dep' : [ref] ,
             'actions' : [ "bwa index %(dependencies)s" ] 
             }
Arg = str # (there is probably a smarter way do this)
# need a special action type to accept the argument hmm
#GAction = Union[str, Callable[[List[FPath], List[FPath]], ActionReturn]]
#ActionReturn = Union[Failure, None, Dict[str, bool]]

def task_bwa_test() -> BoolArgJob:
    #note: unpaired_reads is a special arg
    def bwa_mem_action(targets: List[FPath], dependencies: List[FPath], unpaired_reads: bool) -> None:
        r1, r2, u1, u2 = targets
        paired_bam = Bam('paired.bam')
        if unpaired_reads:
            unpaired_bam = Bam('unpaired.bam')
        else:
            pass 
        subprocess.call(f"touch {targets}")
        return None
            
    return { 'targets' : [unsorted_bam],
             'file_dep' : index_files + [trimmed_r1, trimmed_r2, up_trimmed1, up_trimmed2], 
             'actions' : [ bwa_mem_action ],
             'getargs'   : { 'unpaired_reads' : ('trimmomatic', 'unpaired_reads')} }
             


def task_sort_bam() -> Job:
    # how can these targets and dependencies be parameterizable like they are in, say, `dagr`?
    return { 'targets' : [sorted_bam],
            'file_dep' : [unsorted_bam],
            'actions'  : [ "samtools sort %(dependencies)s -f %(targets)s" ] }

def task_index_bam() -> Job: 

    return { 'targets'  : [index_bai],
             'file_dep' : [sorted_bam],
             'actions'  : [ "samtools index %(dependencies)s" ] }

def task_fqstats() -> Job:
    cfg: Config = yaml.load(open(args['config']))
    sample = 'SAMPLE' # args.get('sample', 'SAMPLE')
    skip_stats = cfg['fqstats']['skip']['default'] # or accept commandline arg 
    if skip_stats:
        return { 'targets' : [], 'file_dep' : [], 'actions': []}
    else: 
        return { 'targets' : [FPath(sample+ '_fqstat').png],
                'file_dep' : [trimmed_r1, trimmed_r2], 
                'actions' : ['fqstats -o %(targets)s %(dependencies)s'] }
# cmd1 = 'lofreq call-parallel --pp-threads ' + str(threads) + ' -f {reference} {bamfile} -o {vcf} ' + lofreq_options

def task_lofreq_index_ref() -> Job:

    return { 'targets' : [ref.fai],
            'actions' : ["lofreq faidx %(dependencies)s"],
            'file_dep' : [ ref ] }
def task_call_variants() -> Job:
 caller = 'lofreq'
 threads = 2
 lofreq_vcf = VCF('lofreq.vcf')

 if caller is 'lofreq':
  return { 'targets' : [lofreq_vcf],
           'file_dep' : [sorted_bam, ref, ref.fai, index_bai],
           'actions' : [ # TODO: must bea  better way to give it individual dependency files
           f"lofreq call-parallel --pp-threads {threads} -f {ref} {sorted_bam} -o %(targets)s" # + lofreq-options
           ] }

 else:
     return cast(Job, {})




