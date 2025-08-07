#!/usr/bin/bash
#SBATCH --time=6:00:00 

SNAKEFILE=$1

if [[ ! -e ${SNAKEFILE} ]]; then
    echo; echo Please add Snakemake file name; echo
    exit
fi    

module load snakemake

source ~/bin/myconda

### Run pipeline
snakemake --profile ./config/snakemake_profile --snakefile ${SNAKEFILE} -p
 

