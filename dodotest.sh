#!/bin/bash
cp ngs_doit/ngs_doit.py dodotest/dodo.py
cd dodotest && doit config=../config.yaml.default r1=R1.fq r2=R2.fq ref=ref.fasta sample=69420 outdir=./
