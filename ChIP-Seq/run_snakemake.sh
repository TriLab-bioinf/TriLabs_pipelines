#!/usr/bin/bash
#SBATCH --time=10-00:00:00 

SNAKEFILE=./workflow/ChIP-Seq_pipeline.smk
PARAMS="$*"

if [[ ! -e ${SNAKEFILE} ]]; then
    echo; echo Please enter path to snakemake file; echo 
    exit
fi    

module load snakemake

source ~/bin/myconda

### Run pipeline
snakemake --profile ../snakemake_profile --snakefile ${SNAKEFILE} -p ${PARAMS}