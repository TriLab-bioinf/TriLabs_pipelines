rule filter_variants:
    input:
        snps = "results/10-subset-snps/all_SNP.vcf.gz",
        indels = "results/11-subset-indels/all_INDELS.vcf.gz",
        ref = f"{genome}"
    output:
        snps = "results/12-filter-snps/all_SNP.filter.vcf.gz",
        indels = "results/13-filter-indels/all_INDELS.filter.vcf.gz"
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
