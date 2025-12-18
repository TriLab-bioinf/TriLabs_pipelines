
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

# Joined gentyping using GenomicsDB 

rule genotype_gvcfs:
    input:
        genomicsdb = rules.genomics_db_import.output.db,
        ref = f"{genome}"
    output:
        vcf = "results/8-all-calls/all.vcf.gz"
    log:
        logs = "logs/8-all-calls/genotypegvcfs.log"
    params:
        exome_intervals = f"--intervals {intervals}" if config["exomeseq"] == True else "",
        exome_interval_padding = f"--interval-padding {config['interval_padding']}" if config["exomeseq"] == True else "",
        extra = "-A StrandBiasBySample",  # optional
        java_opts = "", # optional
    resources:
        mem_mb = 8 * 1024,
        runtime = 24 * 60,
        partition = "normal"
    benchmark:
        "benchmarks/all-calls/genotypegvcfs.tsv"
    shell:
        """
        module load GATK/4.6.0.0

        gatk --java-options "-Xmx4g" GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{input.genomicsdb} \
            -O {output.vcf} \
            {params.extra} {params.exome_intervals} {params.exome_interval_padding} > {log.logs} 2>&1

        """

rule select_variants:
    input:
        vcf = "results/8-all-calls/all.vcf.gz"
    output:
        snps = "results/9-subset-snps/all_SNP.vcf.gz",
        indels = "results/10-subset-indels/all_INDELS.vcf.gz"
    log:
        snps = "logs/subset/snps.log",
        indels = "logs/subset/indels.log"
    params:
        extra = ""
    resources:
        mem_mb = 1024 * 32,
        runtime = 2 * 60,
        partition = "quick"
    benchmark:
        "benchmarks/subset-variants/subset.tsv"
    shell:
        """
        module load GATK/4.6.0.0

        gatk SelectVariants \
            -V {input.vcf} \
            -select-type SNP \
            -O {output.snps}  > {log.snps} 2>&1

        gatk SelectVariants \
            -V {input.vcf} \
            -select-type INDEL \
            -O {output.indels}  > {log.indels} 2>&1
        """

rule filter_variants:
    input:
        snps = "results/9-subset-snps/all_SNP.vcf.gz",
        indels = "results/10-subset-indels/all_INDELS.vcf.gz",
        ref = f"{genome}"
    output:
        snps = "results/11-filter-snps/all_SNP.filter.vcf.gz",
        indels = "results/12-filter-indels/all_INDELS.filter.vcf.gz"
    log:
        snps = "logs/filter_variant/snps.log",
        indels = "logs/filter_variant/indels.log"
    params:
        extra = "",
        java_options = "-Xms4G -Xmx4G -XX:ParallelGCThreads=2"
    resources:
        mem_mb = 1024 * 32,
        runtime = 2 * 60,
        partition = "quick"
    benchmark:
        "benchmarks/subset-variants/subset.tsv"
    shell:
        """
        module load GATK/4.6.0.0

        gatk --java-options '{params.java_options}' VariantFiltration \
            --R {input.ref} \
            --V {input.snps} \
            --window 35 \
            --cluster 3 \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O {output.snps} > {log.snps} 2>&1

        gatk --java-options '{params.java_options}' VariantFiltration \
            --R {input.ref} \
            --V {input.indels} \
            --window 35 \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "FS > 200.0" --filter-name "FS200" \
            -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
            -O {output.indels} > {log.indels} 2>&1
        """
