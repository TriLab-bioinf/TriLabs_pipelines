
if config["count_duplicates"] == True:
    remove_flag = ""
else:
    remove_flag = "--ignoreDup"

# Count reads per feature
rule read_counts:
    input: 
        genome = f"{genome}",
        gtf = f"{annotation}",
        bam = lambda wildcards: expand("results/3-noduplicates/{s}.dedup.bam", s=samples) if (config["flag_dup"] == True) else expand("results/2-map_reads/{s}.Aligned.sortedByCoord.out.bam", s=samples)
    output: counts = "results/4-read_counts/read_counts",
            summary = "results/4-read_counts/read_counts.summary"
    params: f"{config['feat_counts_param']} {remove_flag}"
    benchmark:
        "benchmarks/4-read_counts/counts.tsv"
    threads: 16
    resources:
        cpus_per_task = 16,
        partition = "norm",
        time = "14:00:00",
        mem_mb = 32000
    log: 
        logfile = "logs/4-read_counts/featureCounts.log"
    shell:
        """
        module load subread/2.0.6

        featureCounts {params} -G {input.genome} -T {threads}\
         -a {input.gtf} \
         -o {output.counts} {input.bam} > {log.logfile} 2>&1
        """