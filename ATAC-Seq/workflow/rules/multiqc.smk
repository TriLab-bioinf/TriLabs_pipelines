# Run multiqc report
rule multiqc:
    input: 
        i1 = expand("results/1-trim/{sample}.log", sample=metadata['sample']),
        i2 = expand("results/2-map_reads/{sample}.bam", sample=metadata['sample']),
        i3 = lambda wildcards: expand("results/3-dedup/{sample}.dedup.bam", sample=metadata['sample']) if (config["flag_dup"] == True) else [],
        i4 = expand("results/4-peak_calling/{sample}_peaks.narrowPeak", sample=metadata['sample']),
        i5 = expand("results/4-peak_calling/{sample}_treat_pileup.bw", sample=metadata['sample']),
        i6 = expand("results/4-peak_calling/{sample}_control_lambda.bw", sample=metadata['sample']),
        i7 = expand("results/5-bigwig/{sample}.bw", sample=metadata['sample'])
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