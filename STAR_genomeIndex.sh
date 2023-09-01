#!/bin/bash
#SBATCH --job-name=STAR_genomeIndex --output=%x.%j.out
#SBATCH --mail-type=END,FAIL --mail-user=sls366@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=100G

#-----------------------------------------------------------------------------#
# This script makes index of the ref genome using STAR #
#-----------------------------------------------------------------------------#

module load star/2.7.1a

#- Set variables ----------------------------------------------------------------#

refgen_dir=/home/sls366/prapae/genome_files/
refgen_index=/home/sls366/prapae/genome_files/STAR_index/


#- RUN STAR----------------------------------------------------------------#

STAR --runMode genomeGenerate \
        --genomeDir ${refgen_index} \
        --genomeFastaFiles ${refgen_dir}/GCF_905147795.1_ilPieRapa1.1_genomic.fna \
        --sjdbGTFfile ${refgen_dir}/genomic.gff

#- FIN -----------------------------------------------------------------------#

