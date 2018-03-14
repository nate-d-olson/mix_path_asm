#!/usr/bin/sh

QUAST=/scratch/nolson/quast-4.4/quast.py
for i in {1..17};
  do
    cd NIST_${i}-spades
    /scratch/nolson/quast-4.4/quast.py scaffolds.fasta
    cd ../
  done
