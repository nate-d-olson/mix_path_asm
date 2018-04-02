## Mixed Pathogen RM Assemblies
Scripts used to generate draft assemblies for NIST Mix Pathogen RM

### Initial Spades Assembly
Canu requires genome size estimates, since the reference genomes are not provided/ known.
Will generate a spades assembly to get genome size estimate.

The Spades assembly recommend initial trimming and quality filtering.
Skipping qc step as this assembly is only used for estimating genome size.

__Spades assembly script__

```
#!/usr/bin/sh

SPADES=/scratch/nolson/SPAdes-3.11.1-Linux/bin/spades.py
READDIR=../IGS-PacBio-Data
for i in {1..17};
  do
    fq1=${READDIR}/NISTRM_${i}/ILLUMINA_DATA/*_R1.fastq.gz
    fq2=${READDIR}/NISTRM_${i}/ILLUMINA_DATA/*_R2.fastq.gz
    $SPADES --meta -m 100 -1 ${fq1} -2 ${fq2} -o NIST_${i}-spades
  done
```

Assembly statistics calculated using quast.
run from within Spades assembly directories
```
#!/usr/bin/sh

QUAST=/scratch/nolson/quast-4.4/quast.py
for i in {1..17};
  do
    cd NIST_${i}-spades
    /scratch/nolson/quast-4.4/quast.py scaffolds.fasta
    cd ../
  done
```

QUAST was run on Spade assemblies, `nohup bash run_quast_spades.sh > quast_spades_asm.log &`

## Canu Assembly

Assembly sizes based on Spades results.
Canu was run using `nohup bash run_canu.sh > canu_asm.log &`.

```
#!/usr/bin/sh

CANU=canu-1.7/Linux-amd64/bin/canu
for i in {1..17};
  do
    REPORT=NIST_${i}-spades/quast_results/latest/report.tsv
    GS=$(awk 'BEGIN {FS="\t"}; /Total length\t/ {print $2}' ${REPORT})
    echo ${GS}
    FQ=../IGS-PacBio-Data/NISTRM_${i}/PACBIO_DATA/*/Analysis_Results/Filtered/*fastq
    ${CANU}  -assemble \
        Threads=30 \
        Memory=80 \
        -d NIST_${i}-pacbio \
        -p canu \
        genomeSize=${GS} \
        -pacbio-raw

  done
```

## Pilon Commands
Pilon run commands based on tutorial from http://sepsis-omics.github.io/tutorials/modules/cmdline_assembly/
Running Pilon on assemblies using the following script.

```
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

export -f pilon_pipe
parallel -j 4 pilon_pipe ::: {5..17}
```
