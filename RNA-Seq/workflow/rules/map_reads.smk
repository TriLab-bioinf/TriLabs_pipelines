
# Map reads to reference with STAR

# Function to calculate genome length from a fasta file
import os
from math import log2
def calculate_genome_length(fasta_path):
    length = 0
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Fasta file not found: {fasta_path}")
    with open(fasta_path, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                length += len(line.strip())
    return length

genomeSAindexNbases = min(log2(calculate_genome_length(genome)) / 2, 14)

# Set STARdb path and check if STARdb already exists, otherwise create one
if (stardb_path == "None") and (not os.path.exists("data/stardb/SA")):
    db_path: str = "data/stardb/"

    # Make STARdb
    rule make_star_db:
        input: gen = f"{genome}", 
            ann = f"{annotation}" 
        output: db = "data/stardb/SA"
        resources: 
            mem_mb = 1024 * 64,
            partition = "norm",
            runtime = 48 * 60,
            disk_mb= 1024 * 20
        threads: 16
        params: so = f"{stardb_overhang}",
                db_dir = "data/stardb",
                small_genome = f"--genomeSAindexNbases {genomeSAindexNbases}" # FOR BACTERIAL (SMALL) GENOMES
        shell:
            """
            module load STAR/2.7.11b

            STAR --runMode genomeGenerate \
            --genomeDir {params.db_dir} \
            --genomeFastaFiles {input.gen} \
            --sjdbGTFfile {input.ann} \
            --runThreadN {threads} \
            --sjdbOverhang {params.so} {params.small_genome}

            touch {output.db}  
            """
elif (stardb_path == "None") and (os.path.exists("data/stardb/SA")):
    # STARdb was already created => skip make_star_db
    db_path: str = "data/stardb/"
else:
    # Use an existent STARdb located somewhere else
    db_path: str = stardb_path

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
        db = lambda wildcards: rules.make_star_db.output.db if ((stardb_path == "None") and (not os.path.exists("data/stardb/SA"))) else []
    output: 
        bam = "results/2-map_reads/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "results/2-map_reads/{sample}.Aligned.sortedByCoord.out.bam.bai"
    threads: 8
    resources:
        partition = "norm",
        runtime = 14 * 60,
        mem_mb = 1024 * 64,
        disk_mb= 1024 * 20
    benchmark:
        "benchmarks/2-map_reads/{sample}.bwa.tsv"
    params: 
        stardb_dir = f"{db_path}",
        prefix = "results/2-map_reads/{sample}.",
        spliced_alignment = "--alignIntronMax 1" if config["spliced_alignment"] == False else ""
    log: 
        logfile = "logs/2-map_reads/{sample}.star.log"
    shell:
        """
        module load STAR/2.7.11b samtools/1.21
            
        STAR --runMode alignReads \
                --runThreadN {threads} \
                --genomeDir {params.stardb_dir} \
                --alignSJDBoverhangMin 1 \
                --alignSJoverhangMin 5 \
                --outFilterMismatchNmax 2 \
                --alignEndsType EndToEnd \
                --readFilesIn {input.fq1} {input.fq2} \
                --readFilesCommand zcat \
                --outFileNamePrefix {params.prefix} {params.spliced_alignment} \
                --quantMode GeneCounts \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattrRGline ID:$$ SM:{wildcards.sample} PL:ILLUMINA \
                --outSAMattributes All > {log.logfile} 2>&1

        samtools index -@ 8 {output.bam} >> {log.logfile} 2>&1
        """
