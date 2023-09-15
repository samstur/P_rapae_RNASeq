#!/bin/bash
#SBATCH --job-name=htseq_count --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sls366@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=96:00:00
#SBATCH --mem=20G

#-----------------------------------------------------------------------------#
# This script uses htseq to count and assign reads to genes #
#-----------------------------------------------------------------------------#

source activate conda_env

#- Set variables ----------------------------------------------------------------#

bam_dir=/home/sls366/prapae/bam_dir
count_dir=/home/sls366/prapae/counts_dir
htseq=/home/sls366/.conda/envs/conda_env/bin/htseq-count
ref=/home/sls366/prapae/genome_files/genomic_notRNA_norRNA.gff # I had to subset the gff to remove all tRNA and rRNA to avoid an error from HTSeq.

#- RUN htseq ----------------------------------------------------------------#

files=(${bam_dir}/*_Aligned.out.bam)
for file in ${files[@]}
do
base=`basename ${file} _Aligned.out.bam`
${htseq} -f bam -r pos -s reverse -t exon -i gene ${bam_dir}/${base}_Aligned.out.bam ${ref} > ${count_dir}/${base}_htseqCount

done

## -s reverse indicates that our sequencing is reverse stranded

#- FIN -----------------------------------------------------------------------#
