# Map reads to reference with bowtie2

# Set bowtiedb path and check if bowtiedb already exists, otherwise create one
if (bowtiedb_path == "None") and (not os.path.exists("data/bowtiedb/genome")):
    db_path: str = "data/bowtiedb/"

    # Make bowtiedb
    rule make_bowtie_index:
        input: gen = f"{genome}"
        output: db = "data/bowtiedb/genome",
        resources: 
            mem_mb = 32000,
            partition = "quick",
            runtime = 14 * 60
        threads: 8
        params: 
        shell:
            """
            module load bowtie

            bowtie2-build --threads {threads} {input.gen} {output.db} 
                
            touch {output.db}
            """

elif (bowtiedb_path == "None") and (os.path.exists("data/bowtiedb/genome")):
    # bowtiedb was already created => skip make_bowtie_index
    db_path: str = "data/bowtiedb/"
else:
    # Use an existent bowtiedb located somewhere else
    db_path: str = bowtiedb_path

# Mapping reads to reference
rule map_reads:
    input:
        fq1 = lambda wildcards: [f"results/1-trim/{wildcards.sample}_{wildcards.type}.P.R1.fastq.gz" if (config["trimming"] == True) else f"{reads_dir}/{get_fastq(wildcards.sample, wildcards.type,'fastq_1')}"],
        fq2 = lambda wildcards: (
                [
                    f"results/1-trim/{wildcards.sample}_{wildcards.type}.P.R2.fastq.gz" if (config["trimming"] == True) else f"{reads_dir}/{get_fastq(wildcards.sample, wildcards.type,'fastq_2')}"
                ]
                if config["paired"] == True
                else []
            ),
        db = lambda wildcards: rules.make_bowtie_index.output.db if ((bowtiedb_path == "None") and (not os.path.exists("data/bowtiedb/genome"))) else []
    output: 
        bam = "results/2-map_reads/{sample}_{type}.bam",
        bai = "results/2-map_reads/{sample}_{type}.bam.bai"
    threads: 16
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 1024 * 64,
        disk_mb= 1024 * 20
    benchmark:
        "benchmarks/2-map_reads/{sample}_{type}.bowtie2.tsv"
    params: 
        prefix = "{sample}_{type}",
        bowtiedb_dir = f"{db_path}",
        db_base = "genome"
    log: "logs/2-map_reads/{sample}_{type}.bowtie2.log"        
    shell:
        """
        module load bowtie samtools

        bowtie2 -x {params.bowtiedb_dir}/{params.db_base} \
        --threads {threads} \
        --phred33 \
        --rg-id {params.prefix} --rg SM:{params.prefix} --rg LB:library --rg PL:ILLUMINA \
        -U {input.fq1} \
        --sensitive | samtools sort -@ 8 -O BAM -o {output.bam} -
        
        samtools index -@ 8 {output.bam}
        """