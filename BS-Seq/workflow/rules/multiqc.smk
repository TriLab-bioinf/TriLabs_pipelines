# Run multiqc report

rule multiqc:
    input: 
        i1 = lambda wildcards: expand("results/1-trim/{s}.P.R{r}.fastq.gz", s = samples, r = [1,2]),
        i2 = expand("results/2-bismark/{s}/bismark.bam", s = samples),
        i3 = expand("results/3-noduplicates/{s}.deduplicated.bam", s = samples),
        i4 = expand("results/4-meth-extractor/{s}.deduplicated_splitting_report.txt", s = samples),
        i5 = expand("results/5-bismark-reports/{s}/bismark_report.html", s = samples)
        #i5 = rules.bismark2report.output.bm_report_dir
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