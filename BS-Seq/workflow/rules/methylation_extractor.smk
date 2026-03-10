# Flag duplicated reads in BAM1

rule meth_extractor:
    input:
        dedup_bam = "results/3-noduplicates/{sample}.deduplicated.bam",
        genome = f"{config['bismark']['path']}"
    output:
        split = "results/4-meth-extractor/{sample}.deduplicated_splitting_report.txt",
        mbias_report = "results/4-meth-extractor/{sample}.deduplicated.M-bias.txt"
    threads: 30
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 64000,
        disk_mb= 20480
    params: # place holder for other parameters to pass to the methylation extraction command
        other_meth_params = f"{config['bismark']['other_meth_extract_params']}",
        outdir = "results/4-meth-extractor/",
        processes = 10 # Should be threads/3
    log: 
        logfile = "logs/4-meth-extractor/{sample}.meth-extractor.log"
    benchmark:
        "benchmarks/4-meth-extractor/{sample}.tsv"
    shell:
        """
        module load bismark/0.25.0 samtools/1.15.1

        bismark_methylation_extractor -p --no_overlap \
            --comprehensive --gzip --bedGraph \
            --cytosine_report \
            --genome_folder {input.genome} \
            --multicore {params.processes} \
            --output_dir {params.outdir} \
            {input.dedup_bam} > {log.logfile} 2>&1
        """

rule bismark2report:
    input:
        bm_report = "results/2-bismark/{sample}/bismark_report.txt",
        split_report = "results/4-meth-extractor/{sample}.deduplicated_splitting_report.txt",
        dedup_report = "results/3-noduplicates/{sample}.deduplication_report.txt",
        mbias_report = "results/4-meth-extractor/{sample}.deduplicated.M-bias.txt"
    output: 
        bm_report_dir = directory("results/5-bismark-reports/{sample}"),
        bm_report_html = "results/5-bismark-reports/{sample}/bismark_report.html"
    params: # place holder for other parameters to pass to the methylation extraction command
        other_b2r_params = f"{config['bismark']['other_b2r_params']}"
    threads: 8
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 64000,
        disk_mb= 20480
    log: 
        logfile = "logs/5-bismark-reports/{sample}.bismark-report.log"
    benchmark:
        "benchmarks/5-bismark-reports/{sample}.tsv"
    shell:
        """
        module load bismark/0.25.0 samtools/1.15.1

        # Prepare reports
        bismark2report --dir {output.bm_report_dir} \
            --alignment_report {input.bm_report} \
            --dedup_report {input.dedup_report} \
            --splitting_report {input.split_report}\
            --mbias_report {input.mbias_report} > {log.logfile} 2>&1
        """