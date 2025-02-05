#!/bin/bash
#SBATCH --job-name=upload --output=z.upload
#SBATCH --mail-type=END,FAIL --mail-user=sls366@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=10:00:00
#SBATCH --mem=10GB
 
# Transfer tar file from the google bucket to my HPC directory
gsutil cp gs://gu-biology-pi-paa9/sls366/p_rapae_RNASeq/XPAGU.20230731_210815.ILLUMINA_DATA.1-of-1.tar /home/sls366/prapae/raw_dir/new_raw_dir/
