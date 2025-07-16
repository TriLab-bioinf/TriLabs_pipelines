#!/usr/bin/bash

set -o errexit 

Help()
{
   # Display Help
   echo "Script to print out final shell commands run by the snakemake pipeline. It is necessary to specify the path to the main snakemake file and the shell tools to be output separated by commas."
   echo
   echo "Syntax: ${0} [-p|-c|-o|-h]"
   echo "options:"
   echo "-p = path to snakemake file"
   echo "-c = comma-separated list of the commands to be output (not case sensitive. e.g. module,STAR,samtools,bamCoverage,featureCounts,multiqc )"
   echo "-o = output file [default = ./workflow/snakemake_shell_commands.txt]"
   echo "-h = Prints this help."
   echo
}

# Deals with NULL entry
if [[ -z $1 ]]; then
    Help
    exit 1  
fi    

# Initialize variables with default values
OUTPUT=workflow/snakemake_shell_commands.txt

# Process options
while getopts "hp:c:o:" option; do
   case $option in
        h) # display Help
                Help
                exit;;
        p) SMKFILE=${OPTARG};;
        c) CMD_LIST=${OPTARG};;
        o) OUTPUT=${OPTARG};;
       \?) # incorrect option
           echo "Error, Invalid option"
           echo
             Help
           exit 1
           ;;
       :) echo "option -$OPTARG requires an argument." 
            Help
         exit 1
         ;;   
   esac
done

if [[ ! -e ${SMKFILE} ]]; then
    echo; echo Error, ${SMKFILE} snakemake file not found
fi

OUTPUTDIR=$(dirname ${OUTPUT})
mkdir -p ${OUTPUTDIR}

CMDS=$(echo ${CMD_LIST} | sed 's/,/\\\|/g')

module -q load snakemake

snakemake --snakefile ${SMKFILE} -p -n --forceall | \
    grep -i " (${CMDS})" | \
    grep -v ': \|^rule\|all,' | \
    sed 's/^ \+//' | \
    sed 's/$/\n----------------------/'| \
    sed -e 's/ -/\n  -/g' > ${OUTPUT}

