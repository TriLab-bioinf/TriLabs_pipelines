# Map reads to reference with bwa2

# Set bwadb path and check if bwadb already exists, otherwise create one
if (bwadb_path == "None") and (not os.path.exists("data/bwadb/genome.fa")):
    db_path: str = "data/bwadb/"

    # Make bwadb
    rule make_bwa_index:
        input: gen = f"{genome}"
        output: db = "data/bwadb/genome.fa",
        resources: 
            partition = "norm",
            runtime = 14 * 60,
            mem_mb = 1024 * 100,
            disk_mb= 1024 * 50
        threads: 8
        params: 
        shell:
            """
            module load bwa-mem2/2.2.1

            bwa-mem2 index {input.gen} -p {output.db} 
                
            touch {output.db}
            """

elif (bwadb_path == "None") and (os.path.exists("data/bwadb/genome.fa")):
    # bwadb was already created => skip make_bwa_index
    db_path: str = "data/bwadb/"
else:
    # Use an existent bwadb located somewhere else
    db_path: str = bwadb_path

# Mapping reads to reference
rule map_reads:
    input:
        fq1 = lambda wildcards: [f"results/1-trim/{wildcards.sample}.P.R1.fastq.gz" if (config["trimming"] == True) else f"data/reads/{fq_1[wildcards.sample]}"],
        fq2 = lambda wildcards: (
                [
                    f"results/1-trim/{wildcards.sample}.P.R2.fastq.gz" if (config["trimming"] == True) else f"data/reads/{fq_2[wildcards.sample]}"
                ]
                if config["paired"] == True
                else []
            ),
        db = lambda wildcards: rules.make_bwa_index.output.db if ((bwadb_path == "None") and (not os.path.exists("data/bwadb/genome.fa"))) else []
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
        "benchmarks/2-map_reads/{sample}.bwa.tsv"
    params: 
        prefix = "{sample}",
        bwadb_dir = f"{db_path}",
        db_base = "genome.fa",
        reads = lambda wildcards, input: f"{input.fq1} {input.fq2}" if config["paired"] else f"{input.fq1}",
        rg = "@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
    log: "logs/2-map_reads/{sample}.bwa.log"        
    shell:
        """
        module load bwa-mem2/2.2.1 samtools/1.21

        bwa-mem2 mem {params.bwadb_dir}/{params.db_base} \
        {params.reads} \
        -t {threads} \
        -R '{params.rg}' \
        | samtools sort -@ {threads} -O BAM -o {output.bam} -
        
        samtools index -@ {threads} {output.bam}
        """