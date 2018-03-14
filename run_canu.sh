#!/usr/bin/sh

CANU=/scratch/nolson/canu-1.7/Linux-amd64/bin/canu
for i in {1..17};
  do
    REPORT=NIST_${i}-spades/quast_results/latest/report.tsv
    GS=$(awk 'BEGIN {FS="\t"}; /Total length\t/ {print $2}' ${REPORT})
    FQ=../IGS-PacBio-Data/NISTRM_${i}/PACBIO_DATA/*/Analysis_Results/Filtered/*fastq
    ${CANU}  -d NIST_${i}-canu -p canu genomeSize=${GS} -pacbio-raw ${FQ}
  done
