# type: ignore
from random import randint
import subprocess
from datetime import datetime
# TODO: dependencies is actually a set, making no guaranteesa bout ourder.
def time():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')
strip = '\n' + '='*10 + '\n'
def task_trimmer():
    def run_trim(targets, dependencies):
        n = randint(0, 99999999) 
        print(targets, dependencies) 

        justinput = [dep for dep in dependencies if dep == 'input.special']
        justother = [dep for dep in dependencies if dep == 'other.txt']
        subprocess.getoutput(f'cp {justinput[0]} {targets[0]}')
        subprocess.getoutput(f'echo "{time()}  trim:  {n}" >> {targets[0]}')
        print("hi other.txt! This is task_trimmer\n", open(justother[0]).read())
        print(f"n: {n}")
        print(strip)
        return { 'unpaired_reads' : n }
    return { 'targets'  : ['trim.txt'],
             'actions'  : [run_trim], 
             'file_dep' : ['input.special', 'other.txt'] }
    #TODO: for some reason touching 'dumb.txt' caused other.txt to get remade I think
def task_other():
    return { 'targets' : ['other.txt'], 'actions' : ['echo "this is other"', f'echo "{time()} otherXXX" > %(targets)s'], 'file_dep' : ['input.special'] }


#TODO: pydoit attempted to build even though the dependencies (file deps) weren't available!
# i.e. when I miss-named task_dumb to dumb_task lul

def task_dumb():
    return { 'targets' : ['dumb.txt'], 'actions' : ['echo "this is dumb task!"', 'cp %(dependencies)s %(targets)s'], 'file_dep' : ['input.special'] }

from doit import create_after
@create_after(executed='trimmer', target_regex='.*\.bam')
def task_vcf():
    nonempty = lambda _: False
    print('smart bwa executing', strip) 
    unpaired = nonempty('up1.fq') or nonempty('up2.fq')
    base_dict = { 'targets' : ['lofreq.vcf'], 
            'actions' : [ 'lofreq call %(dependencies)s  > %(targets)s' ],
            'file_dep' : [ref.fai] }
    if unpaired:
        base_dict['file_dep'] += ['merged.bam']
    else:
        base_dict['file_dep'] += ['paired.bam']
    return base_dict
        
def task_compile_fastq():
    return {'targets' : unpaired_compiled_fastq, 'dependencies' : [up1_fq, up2_fq], 'actions' : [ 'cat %(dependencies)s > %(targets)s' ] }
def task_merge_bam():
    return { 'targets' : ['merged.bam'], 'dependencies' : [paired_bam, unpaired_bam], 'actions' : [ 'samtools merge %(targets)s %(dependencies)s' ] }
#@create_after(executed='trimmer', target_regex='.*\.bam')
def task_bwa_test():
    print('smart bwa task closure', strip)
    yield { 'targets' : [ 'paired.bam' ],
            'dependencies' : [ ref, p1_fq, p2_fq ],
            'actions' : [ 'bwa mem %(dependencies)s > %(targets)s' ] } 
    yield { 'targets' : [ 'unpaired.bam' ], 
            'dependencies' : [ ref, unpaired_compiled_fastq ], 
            'actions' : [ 'bwa mem %(dependencies)s > %(targets)s' ] }
    def smart_bwa(targets, dependencies, unpaired_reads):
        nonempty = lambda _: False
        print('smart bwa executing', strip)
        subprocess.getoutput(f"touch {' '.join(targets)}")
        subprocess.getoutput(f'echo "{time()} Smart:  up: {unpaired_reads}" >> {targets[0]}')
        subprocess.getoutput(f'cat {dependencies[0]} {dependencies[1]} > {targets[1]}')

    return { 'targets'   : [ 'smart-bwa.txt', 'bwa-compile.txt'],
             'file_dep'  : [ 'dumb.txt', 'trim.txt' ],
             'actions'   : [ smart_bwa ],
             'getargs'   :  {'unpaired_reads' : ('trimmer', 'unpaired_reads') }}
