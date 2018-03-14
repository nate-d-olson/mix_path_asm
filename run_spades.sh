#!/usr/bin/bash

SPADES=/scratch/nolson/SPAdes-3.11.1-Linux/bin/spades.py
READDIR=../IGS-PacBio-Data
for i in {1..17};
  do
    fq1=${READDIR}/NISTRM_${i}/ILLUMINA_DATA/*_R1.fastq.gz
    fq2=${READDIR}/NISTRM_${i}/ILLUMINA_DATA/*_R2.fastq.gz
    $SPADES --meta -t 30 -m 100 -1 ${fq1} -2 ${fq2} -o NIST_${i}-spades
  done
