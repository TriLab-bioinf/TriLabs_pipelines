#!/usr/bin/bash
#SBATCH --time=10-00:00:00 

SNAKEFILE=./workflow/GATK_pipeline.smk
PARAMS=$1

if [[ ! -e ${SNAKEFILE} ]]; then
    echo; echo Please enter path to snakemake file; echo 
    exit
fi    

module load snakemake

source ~/bin/myconda

### Run pipeline
snakemake --profile ../snakemake_profile --snakefile ${SNAKEFILE} -p ${PARAMS}
 
