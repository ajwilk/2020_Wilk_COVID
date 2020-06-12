#!/bin/bash
#
#SBATCH --job-name=dropest-hg
#
#SBATCH --time=500:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G

#singularity container for dropEst built from https://hub.docker.com/r/stratust/dropest
#STAR v2.5.4b uaed for this work

export R1="insert_name_of_R1_fastq.fastq.gz"
export R2="insert_name_of_R2_fastq.fastq.gz"
singularity exec /oak/stanford/groups/cblish/aaron/src/stratust-dropest.simg droptag -c /oak/stanford/groups/cblish/aaron/seqwell/config.xml $R1 $R2
ml biology
ml star
ml samtools
export ALIGN=$(echo $R2.* | tr ' ' ,)
STAR --genomeDir=/oak/stanford/groups/cblish/aaron/src/hg37_sars --limitOutSJcollapsed 10000000 --limitIObufferSize=1500000000 --outFilterMultimapNmax=10 --readFilesCommand zcat --readFilesIn $ALIGN
samtools view -S -b Aligned.out.sam > sample.bam
export GTF=/oak/stanford/groups/cblish/aaron/src/hg37_sars/h37.gtf
singularity exec /oak/stanford/groups/cblish/aaron/src/stratust-dropest.simg dropest -V -c /oak/stanford/groups/cblish/aaron/seqwell/config.xml -g $GTF sample.bam
