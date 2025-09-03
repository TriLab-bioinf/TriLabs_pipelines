# Run multiqc report
rule multiqc:
    input: 
        i1 = lambda wildcards: expand("results/1-trim/{s}_{t}.P.R{r}.fastq.gz", s = metadata['sample_ID'], t= metadata['type'], r = [1,2]) if (config["trimming"] == True) else [],
        i2 = expand("results/2-map_reads/{sample}_{type}.bam", sample=metadata['sample_ID'], type=metadata['type']),
        i3 = lambda wildcards: expand("results/3-dedup/{sample}_{type}.dedup.bam", sample=metadata['sample_ID'], type=metadata['type']) if (config["flag_dup"] == True) else [],
        i4 = expand("results/4-peak_calling/{sample}_peaks.narrowPeak", sample=metadata['sample_ID']),
        i5 = expand("results/4-peak_calling/{sample}_treat_pileup.bw", sample=metadata['sample_ID']),
        i6 = expand("results/4-peak_calling/{sample}_control_lambda.bw", sample=metadata['sample_ID']),
        i7 = expand("results/5-bigwig/{sample}_{type}.bw", sample=metadata['sample_ID'], type=metadata['type'])
    output: "results/6-multiqc/multiqc_report.html"
    resources:
        partition = "quick",
        mem_mb = 4000
    benchmark:
        "benchmarks/6-multiqc/multiqc.tsv"
    log: 
        logfile = "logs/6-multiqc/multiqc.log"
    shell:
        """
        module load multiqc
        
        multiqc -f -d -o results/6-multiqc results/ > {log.logfile} 2>&1
        """