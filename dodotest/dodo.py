# type: ignore
from doit import create_after, get_var
import os
import re
from ngs_doit.tools import FqPair
from typing import Optional, List, Tuple, Union, Any, Dict, Callable
from typing_extensions import TypedDict, Literal 
from mypy_extensions import (Arg, DefaultArg, NamedArg,
                             DefaultNamedArg, VarArg, KwArg)
from pathlib import PurePosixPath
from ngs_doit import plotting
from ngs_doit import tools
from ngs_doit.custom_types import (Args,
                                   Job, IdxBamJob, TaskDepJob,
                                   Targets, FileDeps, Actions, 
                                   Failure, PyAction, PathLike, Args)

# TODO: lofreq is python 3.6 so it won't install properly,it's also not in the conda requirements
# also lofreq comes with a suggested snakemake pipline
args: Args = {"config" : get_var('config', ''),
              "outdir" : get_var('outdir', ''),
              "ref"    : get_var('ref',    ''),
              "r1"     : get_var('r1',     ''),
              "r2"     : get_var('r2',     ''),
              "sample" : get_var('sample', '')}

DOIT_CONFIG = {'action_string_formatting': 'old'}

ref = args['ref']
inputR1_fq, inputR2_fq = args['r1'], args['r2']
# pipeline output
filtered1_fq, filtered2_fq = 'filtered.r1.fq', 'filtered.r2.fq'
ref_fai = 'ref.fasta.fai'
up1_fq, up2_fq, p1_fq, p2_fq = ('up1.fq', 'up2.fq', 'p1.fq', 'p2.fq')
index_files =  [ref + suf for suf in  { '.amb', '.sa', '.ann', '.bwt', '.pac'} ]
paired_sam, unpaired_sam, merged_sam, sorted_bam = 'paired.sam', 'unpaired.sam', 'merged.sam', 'sorted.bam'
index_bai = sorted_bam + '.bai'
unpaired_compiled_fastq = 'unpcomp.fastq'
lofreq_vcf = 'lofreq.vcf'
trimmed_r1, trimmed_r2 = 'trimmedR1.fastq', 'trimmedR2.fastq'




#@create_after(executed='trimmomatic', target_regex='.*\.sam') # type: ignore

def task_compile_fastq() -> Job:
    return {'targets' : [unpaired_compiled_fastq], 
            'file_dep' : [up1_fq, up2_fq], 
            'actions' : [ 'cat %(dependencies)s > %(targets)s' ] }

def task_merge_sam() -> Job:
    return { 'targets' : [merged_sam], 
             'file_dep' : [paired_sam, unpaired_sam], 
             'actions' : [ 'samtools merge %(targets)s %(dependencies)s' ] }

def task_lofreq_index() -> Job:
    return { 'targets' : [ref_fai], 
            'file_dep' : [ref], 
            'actions' : [ 'lofreq faidx %(dependencies)s' ] }
    #TODO: below depends on bwa version. needs to match for picard read dist
SAMHEADER = "@HD\tVN:0.7.17-r1188\tGO:none\tSO:coordinate" # allows for samtools merge 
#SAMHEADER = "@HD\tVN:1.4\tGO:none\tSO:coordinate" # allows for samtools merge 

def task_bwa_unpaired() -> Job:
    return { 'targets' : [ unpaired_sam ], 
            'file_dep' : [ ref, unpaired_compiled_fastq] + index_files ,
            'actions' : [ f'bwa mem -H "{SAMHEADER}" {ref} {unpaired_compiled_fastq} > %(targets)s' ] }

def task_bwa_paired() -> Job:
    return { 'targets' : [ paired_sam ],
            'file_dep' : [ ref, p1_fq, p2_fq] +  index_files,
            'actions' : [ f'bwa mem -H "{SAMHEADER}" {ref} {p1_fq} {p2_fq} > %(targets)s' ] } 

def task_bwa_index() -> Job:
    return { 'targets' : index_files,
             'file_dep' : [ref] ,
             'actions' : [ "bwa index %(dependencies)s" ] 
             }


def task_trimmomatic() -> Job:
    q  = 30 # cfg['trim_reads']['q']['default']
    hc = 10 # cfg['trim_reads']['headcrop']['default']
    return {
            'targets'  : [p1_fq, up1_fq, p2_fq, up2_fq],
            'file_dep' : [filtered1_fq, filtered2_fq],
            'actions'  : [f"trimmomatic PE %(dependencies)s %(targets)s LEADING:{q} TRAILING:{q} HEADCROP:{hc} -phred33"]
            }

def get_index(fn: str) -> Optional[str]:
    ''' returns index path (ie _I1_ or _I2_)  or none.'''
    index = re.sub(r'_R([12])_', r'_I\1_', str(fn))
    return index if os.path.exists(index) else None

def task_ngs_filter() -> Job:
    #TODO: get from input 
    minbq, keepNs = 1, True
    inputI1, inputI2 = get_index(inputR1_fq), get_index(inputR2_fq)
    deps = [inputR1_fq, inputR2_fq]
    if inputI1 and inputI2:
        deps += [inputI1, inputI2]
        idxs: Optional[FqPair] = FqPair(inputI1, inputI2)
    else:
        idxs = None
    return { 'targets' : [filtered1_fq, filtered2_fq],
             'file_dep' : deps, 
             'actions' : [(tools.fastq_filter,  # type: ignore
                 [FqPair(inputR1_fq, inputR2_fq), idxs, 
                     FqPair(filtered1_fq, filtered2_fq), minbq, keepNs] ) ]}

def task_fqstats() -> Job:
    fqstats_png = 'fqstats.png'
    #TODO: implement config stuff
    skip_stats = True #cfg['fqstats']['skip']['default'] # or accept commandline arg 
    if skip_stats:
        return { 'targets' : [], 'file_dep' : [], 'actions': []}

    return { 'targets' : [fqstats_png],
             'file_dep' : [inputR1_fq, inputR2_fq], 
             'actions' : ['fqstats -o %(targets)s %(dependencies)s'] }

@create_after(executed='trimmomatic', target_regex='.\*.bam') # type: ignore
def task_sort_bam() -> TaskDepJob:
    nonempty: Callable[[str], bool] = lambda f: os.path.exists(f) and os.stat(f).st_size > 10
    unpaired = nonempty(up1_fq) or nonempty(up2_fq)
    # annotation necessary b/c of decorator
    base_dict: TaskDepJob = { 'targets' : [sorted_bam],
                              'actions' : [ 'samtools sort %(dependencies)s -o %(targets)s' ],
                              'file_dep' : [],
                              'task_dep' : ['trimmomatic']} # necessary b/c it's a delayed task
    if unpaired:
        base_dict['file_dep'] += [merged_sam]
    else:
        base_dict['file_dep'] += [paired_sam]
    return base_dict

def task_read_qual_dist() -> Job:
    qual_score_dist_txt, qual_score_dist_pdf =  'qual_score_dist.txt', 'qual_score_dist.pdf'
    return { 'targets' : [qual_score_dist_txt, qual_score_dist_pdf],
             'file_dep' : [ref, sorted_bam],
             'actions' : [ f"picard QualityScoreDistribution I={sorted_bam} R={ref} \
                              O={qual_score_dist_txt} CHART={qual_score_dist_pdf}   VALIDATION_STRINGENCY=SILENT" ] }


def task_coverge_plot() -> IdxBamJob:
    #TODO: minbq, alloworphans, overlap etc comes comes from config
    # Could simplify to be like length_dist_plot
    #TODO PNG etc only if ocra or w/e installed
    do_orphans, do_overlaps, minbq = True, False, 0 
    return { 'targets' : [ 'qualdepth.json', 'qualdepth.png', 'qualdepth.html' ],
             'file_dep' : [sorted_bam, index_bai],
             #'task_dep' : ['index_bam'],
             'actions' : [(plotting.plot_coverage, [],  # type: ignore
                 { 'orphans'  : do_orphans, 
                   'overlaps' : do_overlaps,
                   'minbq'    : minbq } )] }

def task_length_dist_plot() -> Job:
    return { 'targets'  :  [ 'readlength.html', 'readlength.json'],
             'file_dep' : [ inputR1_fq, inputR2_fq ],
             'actions' :  [ plotting.plot_length_dist ] } # type: ignore
    #TODO: minbq, alloworphans, overlap etc comes comes from config




def task_index_bam() -> Job:
    return { 'targets'  : [index_bai],
             'file_dep' : [sorted_bam],
             'actions'  : [ "samtools index %(dependencies)s" ] }

def task_variant_caller() -> IdxBamJob:
    d: IdxBamJob = { 'targets' : [lofreq_vcf],
          'file_dep' : [ref, sorted_bam, index_bai],
          'actions' : [],
          #'task_dep' : ['index_bam'] 
          }
    caller = 'lofreq'  # config stuff
    if caller == 'lofreq':
        d['actions'] +=  "lofreq call-parallel --pp-threads {threads} -f {ref} {sorted_bam} -o %(targets)s"
        d['file_dep'] +=  [ ref_fai ]
        return d
    else:
        raise NotImplementedError() 

