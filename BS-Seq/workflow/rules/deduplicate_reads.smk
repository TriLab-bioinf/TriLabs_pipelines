# Flag duplicated reads in BAM1

rule mark_duplicates:
    input: 
        bam = "results/2-bismark/{sample}/bismark.bam"
    output: 
        filtered = "results/2-bismark/{sample}/bismark.nonCG_filtered.bam",
        dedup = "results/3-noduplicates/{sample}.deduplicated.bam",
        report = "results/3-noduplicates/{sample}.deduplication_report.txt"
    params: # place holder for other parameters to pass to the deduplication command
        other_dedup_params = f"{config['bismark']['other_dedup_params']}"
    threads: 8
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 64000,
        disk_mb= 20480
    log: 
        logfile1 = "logs/3-noduplicates/{sample}.filter-non-conversion.log",
        logfile2 = "logs/3-noduplicates/{sample}.dedup.log"
    benchmark:
        "benchmarks/3-noduplicates/{sample}.tsv"
    shell:
        """
        module load bismark/0.25.0

        # Filter out reads that fail the non-conversion check
        filter_non_conversion -p {input.bam} > {log.logfile1} 2>&1

        # Remove duplicated reads 
        deduplicate_bismark --bam --output_dir results/3-noduplicates/ --outfile {wildcards.sample} \
            {params.other_dedup_params} -p {output.filtered} > {log.logfile2} 2>&1
        """