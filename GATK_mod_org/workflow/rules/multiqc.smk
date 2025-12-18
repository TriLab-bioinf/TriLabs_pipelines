# Run multiqc report
rule multiqc:
    input: 
        i1 = lambda wildcards: expand("results/1-trim/{s}.P.R{r}.fastq.gz", s = samples, r = [1,2]) if (config["trimming"] == True) else [],
        i2 = expand("results/2-map_reads/{sample}.bam", sample=samples),
        i3 = expand("results/3-dedup/{sample}.dedup.bam", sample=samples),
        i4 = expand("results/4-recal/{s}.grp", s = samples),
        i5 = expand("results/5-recal/{s}.recal.bam", s = samples),
        i6 = expand("results/6-calls-gvfc/{s}.g.vcf.gz", s = samples),
        i7 = rules.genomics_db_import.output.db,
        i7_bkp = rules.genomics_db_import.output.db_bkp,
        i8 = "results/8-all-calls/all.vcf.gz",
        i10 = "results/9-subset-snps/all_SNP.vcf.gz",
        i11 = "results/10-subset-indels/all_INDELS.vcf.gz",
        i13 = "results/11-filter-snps/all_SNP.filter.vcf.gz",
        i14 = "results/12-filter-indels/all_INDELS.filter.vcf.gz"
    output: "results/13-multiqc/multiqc_report.html"
    resources:
        partition = "quick",
        mem_mb = 4000
    log: 
        logfile = "logs/13-multiqc/multiqc.log"
    shell:
        """
        module load multiqc/1.28
        
        multiqc -f -d -o results/12-multiqc results/ > {log.logfile} 2>&1
        """