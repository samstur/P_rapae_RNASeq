#!/bin/bash
#SBATCH --job-name=strandtest_htseq --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sls366@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

#-----------------------------------------------------------------------------#
# This script uses htseq to count and assign reads to genes #
#-----------------------------------------------------------------------------#

# activates a virtual conda environment that I've called conda_env
source activate conda_env

#- Set variables ----------------------------------------------------------------#

bam_dir=/home/sls366/prapae/bam_dir
count_dir=/home/sls366/prapae/counts_dir  
htseq=/home/sls366/.conda/envs/conda_env/bin/htseq-count
ref=/home/sls366/prapae/genome_files/genomic_notRNA_norRNA.gff # I had to subset the gff to remove all tRNA and rRNA to avoid an error from HTSeq.

#- RUN htseq ----------------------------------------------------------------#

# first line runs htseq with strandedness but forward on one file
${htseq} -f bam -r pos -s yes -t exon -i gene ${bam_dir}/XPAGU_20230717_A00904_IL100301767_N4UD-A01_L004_16F2_Aligned.out.bam ${ref} > ${count_dir}/fr_XPAGU_20230717_A00904_IL100301767_N4UD-A01_L004_16F2_htseqCount

# second line runs htseq with strandedness but reverse on one file
${htseq} -f bam -r pos -s reverse -t exon -i gene ${bam_dir}/XPAGU_20230717_A00904_IL100301767_N4UD-A01_L004_16F2_Aligned.out.bam ${ref} > ${count_dir}/rv_XPAGU_20230717_A00904_IL100301767_N4UD-A01_L004_16F2_htseqCount

# third line runs htseq without strandedness on one file
${htseq} -f bam -r pos -s no -t exon -i gene ${bam_dir}/XPAGU_20230717_A00904_IL100301767_N4UD-A01_L004_16F2_Aligned.out.bam ${ref} > ${count_dir}/un_XPAGU_20230717_A00904_IL100301767_N4UD-A01_L004_16F2_htseqCount

# runs htseq on single file with different strandedness so we can which output has the least number of no features

