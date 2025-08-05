# Run multiqc report
rule multiqc:
    input: 
        i1 = expand("logs/1-trim/{sample}.log", sample=metadata['sample']),
        i2 = expand("results/2-map_reads/{sample}.bam", sample=metadata['sample']),
        i3 = lambda wildcards: expand("results/3-dedup/{sample}.dedup.bam", sample=metadata['sample']) if (config["flag_dup"] == True) else [],
        i4 = expand("results/4-bqsr/{sample}.recal.bam", sample=metadata['sample']),
        i5 = expand("results/5-vcf/{sample}.raw.g.vcf.gz", sample=metadata['sample']),
        i6 = "results/6-merged/all.wgs.combined.vcf.gz",
        i7 = "results/7-filtered_vcf/all.wgs.indel.SNP.recalibrated_99.9.vcf.gz"
    output: "results/8-multiqc/multiqc_report.html"
    resources:
        partition = "quick",
        mem_mb = 4000
    benchmark:
        "benchmarks/8-multiqc/multiqc.tsv"
    log: 
        logfile = "logs/8-multiqc/multiqc.log"
    shell:
        """
        module load multiqc
        
        multiqc -f -d -o results/8-multiqc results/ > {log.logfile} 2>&1
        """