# Flag duplicated reads in BAM1

rule mark_duplicates:
    input: 
        bam = "results/2-map_reads/{sample}_{type}.bam"
    output: 
        dedupbam = "results/3-dedup/{sample}_{type}.dedup.bam",
        metrics = "results/3-dedup/{sample}_{type}.metrics.txt",
        tmpdir = temp(directory("results/3-dedup/{sample}_{type}.tmp"))
    params: 
        "--READ_NAME_REGEX null --REMOVE_DUPLICATES false"
    log: 
        logfile = "results/3-dedup/{sample}_{type}.dedup.log"
    threads: 24
    benchmark:
        "benchmarks/3-dedup/{sample}_{type}.tsv"
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 1024 * 100,
        disk_mb= 1024 * 20,
        gres = "lscratch:80"
    shell:
        """
        module load picard samtools

        java -Xmx80g -Xms80g -XX:ParallelGCThreads=24 -jar $PICARDJARPATH/picard.jar   MarkDuplicates \
         -I {input.bam} \
         -O {output.dedupbam} \
         -M {output.metrics} \
         {params} \
         --TMP_DIR {output.tmpdir} > {log.logfile} 2>&1

        samtools index {output.dedupbam}
        """