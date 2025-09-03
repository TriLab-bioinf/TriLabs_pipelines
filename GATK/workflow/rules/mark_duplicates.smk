# Flag duplicated reads in BAM1

rule mark_duplicates:
    input: 
        bam = "results/2-map_reads/{sample}.bam"
    output: 
        dedupbam = "results/3-dedup/{sample}.dedup.bam",
        metrics = "results/3-dedup/{sample}.metrics.txt",
        tmpdir = temp(directory("results/3-dedup/{sample}.tmp"))
    params: 
        "--CREATE_INDEX true"
    log: 
        logfile = "logs/3-dedup/{sample}.dedup.log"
    threads: 24
    benchmark:
        "benchmarks/3-dedup/{sample}.tsv"
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 1024 * 100,
        disk_mb= 1024 * 20,
        gres = "lscratch:80"
    shell:
        """
        module load GATK/4.6.0.0

        gatk MarkDuplicates \
         -I {input.bam} \
         -O {output.dedupbam} \
         -M {output.metrics} \
         {params} \
         --TMP_DIR {output.tmpdir} > {log.logfile} 2>&1
        """