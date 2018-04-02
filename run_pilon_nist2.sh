#!/usr/bin/sh


pilon_pipe() {
  ASM=canu.contigs.fasta
  ASMPATH=/scratch/jackson/mix_path_asm/NIST_${1}-canu/${ASM}
  READDIR=/scratch/jackson/IGS-PacBio-Data
  fq1=${READDIR}/NISTRM_${1}/ILLUMINA_DATA/*_R1.fastq.gz
  fq2=${READDIR}/NISTRM_${1}/ILLUMINA_DATA/*_R2.fastq.gz

  PILON=/scratch/nolson/pilon-1.22.jar

  ## Create working directory
  PILONDIR=/galaxy/NIST_${1}-pilon
  mkdir ${PILONDIR}
  cd ${PILONDIR}
  cp ${ASMPATH} ${ASM}

  ## Index assembly
  bwa index ${ASM}
  samtools faidx ${ASM}

  ## Map reads to assemblies and convert to sorted bam
  bwa mem -t 10 ${ASM} ${fq1} ${fq2} | samtools sort > aln.bam

  ## Index bam
  samtools index aln.bam

  ## Run pilon  
  java -Xmx25G -jar ${PILON} --genome ${ASM} --frags aln.bam \
    --fix all --mindepth 0.5 --changes --threads 10
}

pilon_pipe 2
