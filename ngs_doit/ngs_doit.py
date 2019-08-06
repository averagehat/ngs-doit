# type: ignore
from datetime import datetime
from doit import create_after, get_var
import os
def time():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')
strip = '\n' + '='*10 + '\n'

DOIT_CONFIG = {'action_string_formatting': 'old'}
ref = 'ref.fasta'
ref_fai = 'ref.fasta.fai'
up1_fq, up2_fq, p1_fq, p2_fq = ('up1.fq', 'up2.fq', 'p1.fq', 'p2.fq')
index_files =  [ref + suf for suf in  { '.amb', '.sa', '.ann', '.bwt', '.pac'} ]
paired_sam, unpaired_sam, merged_sam, sorted_bam = 'u.sam', 'p.sam', 'M.sam', 'sorted.bam'
index_bai = sorted_bam + '.bai'
unpaired_compiled_fastq = 'unpcomp.fastq'
lofreq_vcf = 'lofreq.vcf'

def task_lofreq_index_ref(): 
    return { 'targets' : [ref_fai],
            'actions' : ["lofreq faidx %(dependencies)s"],
            'file_dep' : [ ref ] }

def task_variant_caller():
    d = { 'targets' : [lofreq_vcf],
           'file_dep' : [ref, sorted_bam],
           'actions' : [],
            'task_dep' : ['index_bam'
           ] }
    caller = 'lofreq'  # config stuff
    if caller is 'lofreq':
        d['actions'] +=  ["lofreq call-parallel --pp-threads {threads} -f {ref} {sorted_bam} -o %(targets)s"]
        d['task_dep'] +=  ['lofreq_index']
    else:
        raise NotImplementedError() 

def task_index_bam():
    return { 'targets'  : [index_bai],
             'file_dep' : [sorted_bam],
             'actions'  : [ "samtools index %(dependencies)s" ] }

sorted_bam = 'sorted.bam'

@create_after(executed='trimmer', target_regex='.*\.sam')
def task_sort_bam():
    nonempty = lambda f: os.stat(f).st_size > 0 
    unpaired = nonempty(up1_fq) or nonempty(up2_fq)
    base_dict = { 'targets' : [sorted_bam],
                  'actions' : [ 'samtools sort %(dependencies)s -f %(targets)s' ],
                  'file_dep' : [],
                  'task_dep' : ['trimmomatic']} # necessary b/c it's a delayed task
    if unpaired:
        base_dict['file_dep'] += [merged_sam]
    else:
        base_dict['file_dep'] += [paired_sam]
    return base_dict

def task_compile_fastq():
    return {'targets' : [unpaired_compiled_fastq], 'file_dep' : [up1_fq, up2_fq], 'actions' : [ 'cat %(dependencies)s > %(targets)s' ] }

def task_merge_sam():
    return { 'targets' : [merged_sam], 'file_dep' : [paired_sam, unpaired_sam], 'actions' : [ 'samtools merge %(targets)s %(dependencies)s' ] }

def task_lofreq_index():
    return { 'targets' : [ref_fai], 'file_dep' : [ref], 'actions' : [ 'lofreq faidx %(dependencies)s' ] }

def task_bwa_unpaired():
    return { 'targets' : [ unpaired_sam ], 
            'file_dep' : [ ref, unpaired_compiled_fastq ],   # cat -> bwa mem
            'actions' : [ 'bwa mem %(dependencies)s > %(targets)s' ] }

def task_bwa_paired():
    return { 'targets' : [ paired_sam ],
            'file_dep' : [ ref, p1_fq, p2_fq ],  # cat -> bwa mem
            'actions' : [ 'bwa mem %(dependencies)s > %(targets)s' ] } 
    

def task_bwa_index(): 
    return { 'targets' : index_files,
             'file_dep' : [ref] ,
             'actions' : [ "bwa index %(dependencies)s" ] 
             }

filtered1_fq, filtered2_fq = 'filtered.r1.fq', 'filtered.r2.fq'

def task_trimmomatic():
    q  = 30 # cfg['trim_reads']['q']['default']
    hc = 10 # cfg['trim_reads']['headcrop']['default']
    return {
            'targets'  : [p1_fq, up1_fq, p2_fq, up2_fq],
            'file_dep' : [filtered1_fq, filtered1_fq],
            'actions'  : [f"trimmomatic PE %(dependencies)s %(targets)s LEADING:{q} TRAILING:{q} HEADCROP:{hc}"]
            }

def task_ngs_filter():
    # TODO: some way to make this okay 
    def filter_action(targets, dependencies):
        input1_fq, input2_fq = args['R1'], args['R2'] 
        # index quality min (avg? I don't remember)
        # avg qual min
        # N maximum count
        raise NotImplementedError()
    return { 'targets' : [filtered1_fq, filtered2_fq],
             'dependencies' : [input1_fq, input2_fq],
             'actions' : [filter_action] }

def task_fqstats():
    fqstats_png = 'fqstats.png'
    #TODO: implement config stuff
    skip_stats = True #cfg['fqstats']['skip']['default'] # or accept commandline arg 
    if skip_stats:
        return { 'targets' : [], 'file_dep' : [], 'actions': []}
    else: 
        return { 'targets' : [fqstats_png],
                'file_dep' : [trimmed_r1, trimmed_r2], 
                'actions' : ['fqstats -o %(targets)s %(dependencies)s'] }

def task_read_qual_dist():
    qual_score_dist_txt, qual_score_dist_pdf =  'qual_score_dist.txt', 'qual_score_dist.pdf'
    return { 'targets' : [qual_score_dist_txt, qual_score_dist_pdf],
             'file_dep' : [ref, sorted_bam],
             'actions' : f"picard QualityScoreDistribution I={sorted_bam} R={ref} \
                              O={qual_score_dist_txt} CHART={qual_score_dist_pdf}" }

from pathlib import PurePosixPath
from ngs_doit import plotting
from ngs_doit.plotting import PileupOptions, plot_coverage

def task_coverge_plot():
    #TODO: minbq, alloworphans, overlap etc comes comes from config
    basepath = PurePosixPath('qualdepth')
    formats = ['html', 'json'] # + PNG if ocra or w/e installed
    tgts = list(map(basepath.with_suffix, formats))
    do_orphans, do_overlaps, minbq = True, False, 0
    plotfunc = lambda ds, ts: plotting.plot_coverage(ds[0], ts[0].with_suffix(''),
                      format=formats, orphans=do_orphans, 
                      overlaps=do_overlaps, minbq=minbq)
    return { 'targets' : tgts,
             'file_dep' : [sorted_bam],
             'actions' : [plotfunc] }



from typing import Optional
from typing_extensions import TypedDict
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
