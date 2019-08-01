# type: ignore
from datetime import datetime
from doit import create_after
import os
def time():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')
strip = '\n' + '='*10 + '\n'

up1_fq, up2_fq, p1_fq, p2_fq = ('up1.fq', 'up2.fq', 'p1.fq', 'p2.fq')
paired_sam, unpaired_sam, merged_sam = 'u.sam', 'p.sam', 'M.sam'
unpaired_compiled_fastq = 'unpcomp.fastq'
ref = 'ref.fasta'
ref_fai = 'ref.fasta.fai'


def task_trimmer():
    print('trimmerbody')
    def run_trim(targets, dependencies):
        depstr = open(dependencies[0]).read()
        for x in targets:
            with open(x, 'w') as out:
                if 'up' in x:
                    continue
                out.write(depstr + '\n' + x + '\n')
    return { 'targets'  : [p1_fq, p2_fq, up1_fq, up2_fq],
             'actions'  : [run_trim], 
             'file_dep' : ['input.special'] }

@create_after(executed='trimmer', target_regex='.*\.sam')
def task_vcf():
    print('vcf closure X0')
    nonempty = lambda f: os.stat(f).st_size > 0 
    unpaired = nonempty(up1_fq) or nonempty(up2_fq)
    base_dict = { 'targets' : ['lofreq.vcf'], 
            #'actions' : [ 'lofreq call %(dependencies)s  > %(targets)s' ],
            'actions' : [ 'cat %(dependencies)s  > %(targets)s' ],
            'file_dep' : [ref_fai],
            'task_dep' : ['trimmer']} # necessary b/c it's a delayed task
    if unpaired:
        base_dict['file_dep'] += [merged_sam]
    else:
        base_dict['file_dep'] += [paired_sam]
    return base_dict

Job = None
def task_compile_fastq():
    return {'targets' : [unpaired_compiled_fastq], 'file_dep' : [up1_fq, up2_fq], 'actions' : [ 'cat %(dependencies)s > %(targets)s' ] }

def task_merge_sam():   # cat -> samtools merge and drop `>`
    return { 'targets' : [merged_sam], 'file_dep' : [paired_sam, unpaired_sam], 'actions' : [ 'cat %(dependencies)s > %(targets)s' ] }
#TODO: return { 'targets' : [merged_sam], 'file_dep' : [paired_sam, unpaired_sam], 'actions' : [ 'samtools merge %(targets)s %(dependencies)s' ] }

def task_lofreq_index():
    return { 'targets' : [ref_fai], 'file_dep' : [ref], 'actions' : [ 'cp %(dependencies)s %(targets)s', 'echo "INDEXED!" >> %(targets)s']}

#@create_after(executed='trimmer', target_regex='.*\.sam')

def task_bwa_unpaired():
    return { 'targets' : [ unpaired_sam ], 
            'file_dep' : [ ref, unpaired_compiled_fastq ],   # cat -> bwa mem
            'actions' : [ 'cat %(dependencies)s > %(targets)s' ] }

def task_bwa_paired():
    return { 'targets' : [ paired_sam ],
            'file_dep' : [ ref, p1_fq, p2_fq ],  # cat -> bwa mem
            'actions' : [ 'cat %(dependencies)s > %(targets)s' ] } 
    
index_files =  [ref + suf for suf in  { '.amb', '.sa', '.ann', '.bwt', '.pac'} ]
def task_bwa_index() -> Job:

    return { 'targets' : index_files,
             'file_dep' : [ref] ,
             'actions' : [ "bwa index %(dependencies)s" ] 
             }

