# Flag duplicated reads in BAM1

rule mark_duplicates:
    input: 
        bam = "results/2-map_reads/{sample}.Aligned.sortedByCoord.out.bam"
    output: 
        nodupbam = "results/3-noduplicates/{sample}.dedup.bam",
        metrics = "results/3-noduplicates/{sample}.metrics.txt",
        tmpdir = temp(directory("results/3-noduplicates/{sample}.tmp"))
    params: 
        "--REMOVE_DUPLICATES false --CREATE_INDEX true"
    threads: 8
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 64000,
        disk_mb= 20480
    log: 
        logfile = "logs/3-noduplicates/{sample}.dedup.log"
    benchmark:
        "benchmarks/3-noduplicates/{sample}.tsv"
    shell:
        """
        module load GATK/4.6.0.0 picard/3.3.0 

        gatk --java-options "-Xmx64G -Djava.io.tmpdir={output.tmpdir} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" MarkDuplicates \
                {params} \
                -I {input.bam} \
                -O {output.nodupbam} \
                -M {output.metrics} \
                --READ_NAME_REGEX null \
                --TMP_DIR {output.tmpdir} > {log.logfile} 2>&1
        """