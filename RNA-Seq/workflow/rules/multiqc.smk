# Run multiqc report

rule multiqc:
    input: 
        i1 = lambda wildcards: expand("results/1-trim/{s}.P.R{r}.fastq.gz", s = samples, r = [1,2]) if (config["trimming"] == True) else [],
        i2 = expand("results/2-map_reads/{s}.Aligned.sortedByCoord.out.bam", s = samples),
        i3 = lambda wildcards: expand("results/3-noduplicates/{s}.dedup.bam", s = samples) if (config["flag_dup"] == True) else [],
        i4 = "results/4-read_counts/read_counts",
        i5 = expand("results/5-bigwig/{s}.bw", s = samples)
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