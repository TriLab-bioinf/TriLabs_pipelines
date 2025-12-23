
rule haplotype_caller_gvcf:
    input:
        # single or list of bam files
        bam = "results/5-recal/{sample}.recal.bam",
        ref = f"{genome}"
        # known="dbsnp.vcf"  # optional
    output:
        gvcf = "results/6-calls-gvfc/{sample}.g.vcf.gz"
    log:
        "logs/6-calls-gvfc/{sample}.log",
    params:
        exome_intervals = f"--intervals {intervals}" if config["exomeseq"] == True else "",
        exome_interval_padding = f"--interval-padding {config['interval_padding']}" if config["exomeseq"] == True else "",
        extra = f"-A StrandBiasBySample --standard-min-confidence-threshold-for-calling 30 -dont-use-soft-clipped-bases ",  # optional
        java_opts = "",  # optional
    threads: 4
    resources:
        partition = "normal",
        runtime = 8 * 60,
        mem_mb= 8 * 1024
    benchmark:
        "benchmarks/calls-gvfc/{sample}.g.vcf.tsv"
    shell:
        """
        module load GATK/4.6.0.0

        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -ERC GVCF \
            --native-pair-hmm-threads {threads} \
            {params.extra} {params.exome_intervals} {params.exome_interval_padding} > {log} 2>&1
        """

# build variant DB with GenomicsDBImport wrapper

rule genomics_db_import:
    input:
        expand("results/6-calls-gvfc/{s}.g.vcf.gz", s = samples)
    output:
        db = directory("results/7-db"),
        db_bkp = directory("results/7-db-bkp")
    log:
        "logs/7-db/genomicsdbimport.log",
    params:
        joined_inputs = lambda wildcards, input: " --variant ".join(input),
        intervals = f"{whole_chr_intervals}",  # required: intervals to restrict the analysis to a subset of the genome. Mandatory for WES
        db_action = process_db_action(f"{db_action}"),
        optimize_speed = "--bypass-feature-reader true --genomicsdb-shared-posixfs-optimizations true",
        too_many_contigs = "--merge-contigs-into-num-partitions 25" if config["merge_contigs"] == True else "",  # optional
        java_opts = "-Xmx4g -Xms4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true",  # optional
    threads: 4
    resources:
        mem_mb = lambda wildcards, input: max([input.size_mb * 1.6, 8 * 1024]),
        runtime = 48 * 60,
        partition = "normal",
        disk_mb = 10 * 1024
    benchmark:
        "benchmarks/db/genomics-db-import.tsv"
    shell:
        """
        module load GATK/4.6.0.0

        gatk --java-options '{params.java_opts}' GenomicsDBImport \
            --reader-threads {threads} \
            --variant {params.joined_inputs} \
            --intervals {params.intervals} \
            {params.optimize_speed} {params.too_many_contigs} \
            {params.db_action} {output.db} > {log} 2>&1
        
        mkdir -p {output.db_bkp} >> {log} 2>&1
        
        tar -czf {output.db_bkp}/genomeDB_$(echo $$).tar.gz {output.db} >> {log} 2>&1
        
        touch results/7-db/genomeDB_id.$$
        
        """

# Split merged.exons.bed file into chr-based intervals for running genotype_gvcfs rule by chromosome:
rule split_intervals_by_chr:
    input:
        bed = intervals if config["exomeseq"] == True else whole_chr_intervals
    output:
        chr_bed = expand("results/intervals/{chrom}.bed", chrom = get_chr_ids(whole_chr_intervals)) 
    run:
        import os
        os.makedirs("results/intervals/", exist_ok=True)
        chr_interval = dict()
        with open(input.bed, 'r') as bedfile:
            for line in bedfile:
                chrom = line.strip().split('\t')[0]
                chr_interval.setdefault(chrom, []).append(line)
        
        for chrom in chr_interval.keys():
            with open(f"results/intervals/{chrom}.bed", 'w') as f:
                for line in chr_interval[chrom][0:10]:
                    f.write(line)
        
# Joined gentyping using GenomicsDB 
rule genotype_gvcfs_by_chr:
    input:
        genomicsdb = rules.genomics_db_import.output.db,
        ref = f"{genome}",
        intvl = "results/intervals/{chrom}.bed"
    output:
        vcf = "results/8-all-calls-by-chr/all_samples.{chrom}.vcf.gz"
    log:
        logs = "logs/8-all-calls-by-chr/genotypegvcfs.{chrom}.log"
    params:
        exome_intervals = f"--intervals",
        exome_interval_padding = f"--interval-padding {config['interval_padding']}" if config["exomeseq"] == True else "",
        extra = "-A StrandBiasBySample",  # optional
        java_opts = "", # optional
    resources:
        mem_mb = 8 * 1024,
        runtime = 2 * 24 * 60,
        partition = "normal"
    benchmark:
        "benchmarks/all-calls-by-chr/genotypegvcfs.{chrom}.tsv"
    shell:
        """
        module load GATK/4.6.0.0

        gatk --java-options "-Xmx4g" GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{input.genomicsdb} \
            -O {output.vcf} \
            {params.extra} \
            {params.exome_intervals} {input.intvl} \
            {params.exome_interval_padding} > {log.logs} 2>&1

        """