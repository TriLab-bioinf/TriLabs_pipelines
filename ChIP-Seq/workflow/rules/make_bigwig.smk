if config["bw_ignore_dup"] == True:
    dup_flag = "--ignoreDuplicates"
else:
    dup_flag = ""

# Make bigwig files from bam files 
rule make_bigwig:
    input: lambda wildcards: "results/3-dedup/{sample}_{type}.dedup.bam" if (config["flag_dup"] == True) else "results/2-map_reads/{sample}_{type}.bam"
    output: "results/5-bigwig/{sample}_{type}.bw"
    params: f"--binSize {config['bw_bin']} --normalizeUsing {config['bw_normalization']} {dup_flag}"
    threads: 8
    resources:
        partition = "norm",
        runtime = 8 * 60
    benchmark:
        "benchmarks/5-bigwig/{sample}_{type}.bw.tsv"
    log: 
        logfile = "logs/5-bigwig/{sample}_{type}.bamCoverage.log"
    shell:
        """
        module load deeptools/3.5.6

        bamCoverage -p {threads} -b {input} -o {output} {params} > {log.logfile} 2>&1
        """