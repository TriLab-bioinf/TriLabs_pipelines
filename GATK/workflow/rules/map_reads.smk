# Map reads to reference with bwa2

# Set bwadb path and check if bwadb already exists, otherwise create one
if (bwadb_path == "None") and (not os.path.exists("data/bwadb/genome")):
    db_path: str = "data/bwadb/"

    # Make bwadb
    rule make_bwa_index:
        input: gen = f"{genome}"
        output: db = "data/bwadb/genome",
        resources: 
            mem_mb = 32000,
            partition = "quick",
            runtime = 14 * 60
        threads: 8
        params: 
        shell:
            """
            module load bwa/0.7.17

            bwa index {input.gen} -p {output.db} 
                
            touch {output.db}
            """

elif (bwadb_path == "None") and (os.path.exists("data/bwadb/genome")):
    # bwadb was already created => skip make_bwa_index
    db_path: str = "data/bwadb/"
else:
    # Use an existent bwadb located somewhere else
    db_path: str = bwadb_path

# Mapping reads to reference
rule map_reads:
    input:
        fq1 = "results/1-trim/{sample}_R1.paired.fastq.gz",
        fq2 = "results/1-trim/{sample}_R2.paired.fastq.gz",
        db = lambda wildcards: rules.make_bwa_index.output.db if ((bwadb_path == "None") and (not os.path.exists("data/bwadb/genome"))) else []
    output: 
        bam = "results/2-map_reads/{sample}.bam",
        bai = "results/2-map_reads/{sample}.bam.bai"
    threads: 16
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 1024 * 64,
        disk_mb= 1024 * 20
    benchmark:
        "benchmarks/2-map_reads/{sample}.bwa2.tsv"
    params: 
        prefix = "{sample}",
        bwadb_dir = f"{db_path}",
        db_base = "data/bwadb/genome",
        rg = "@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
    log: "logs/2-map_reads/{sample}.bwa2.log"        
    shell:
        """
        module load bwa/0.7.17 samtools/1.21

        bwa mem {params.db_base} \
        {input.fq1}  {input.fq2} \
        -t {threads} \
        -R '{params.rg}' \
        | samtools sort -@ {threads} -O BAM -o {output.bam} -
        
        samtools index -@ {threads} {output.bam}
        """