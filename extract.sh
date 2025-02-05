#!/bin/bash
#SBATCH --job-name=extract --output=z.extract
#SBATCH --mail-type=END,FAIL --mail-user=sls366@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=10:00:00
#SBATCH --mem=10MB

tar -xf XPAGU.20230731_210815.ILLUMINA_DATA.1-of-1.tar
