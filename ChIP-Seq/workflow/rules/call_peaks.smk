# Call peaks with macs2
rule peak_calling:
    input:
        treat = "results/3-dedup/{sample}_treat.dedup.bam",
        control = "results/3-dedup/{sample}_input.dedup.bam"
    output:
        narrow_peak = "results/4-peak_calling/{sample}_peaks.narrowPeak" if config["peak_calling"]["broad"] == False else [],
        broad_peak = "results/4-peak_calling/{sample}_peaks.broadPeak" if config["peak_calling"]["broad"] else [],
        gapped_peak = "results/4-peak_calling/{sample}_peaks.gappedPeak" if config["peak_calling"]["broad"] else [],
        control_bdg = "results/4-peak_calling/{sample}_control_lambda.bdg",
        chip_bdg = "results/4-peak_calling/{sample}_treat_pileup.bdg",
    params:
        gsize = config["reference"]["gsize"],
        name = "{sample}",
        outdir = "results/4-peak_calling",
        broad = "--broad" if config["peak_calling"]["broad"] else "",
        misc_opts = config["peak_calling"]["other"],
        qvalue = config["peak_calling"]["q_value"]
    log:
        logmacs2 = "logs/4-peak_calling/{sample}.log"
    shell:
        """
        module load macs

        macs2 callpeak -t {input.treat} -c {input.control} \
            -f BAM -g {params.gsize} \
            -n {params.name} \
            -B -q {params.qvalue} \
            --outdir {params.outdir} {params.broad} {params.misc_opts} > {log.logmacs2} 2>&1

        """

# convert bedgraph to bigwig
rule bigwig:
    input: 
        control = "results/4-peak_calling/{sample}_control_lambda.bdg",
        chip = "results/4-peak_calling/{sample}_treat_pileup.bdg"
    output:
        control_bw = "results/4-peak_calling/{sample}_control_lambda.bw",
        chip_bw = "results/4-peak_calling/{sample}_treat_pileup.bw"
    params:
        idx = f"{idx}"
    log:
        logbigwig_input = "logs/4-peak_calling/{sample}_bigwig_input.log",
        logbigwig_chip = "logs/4-peak_calling/{sample}_bigwig_chip.log"
    shell:
        """
        module load ucsc
        
        cut -f1,2 {params.idx} > data/chrom.sizes
        bedGraphToBigWig {input.control} data/chrom.sizes {output.control_bw} >  {log.logbigwig_input} 2>&1 
        bedGraphToBigWig {input.chip} data/chrom.sizes {output.chip_bw} > {log.logbigwig_chip} 2>&1
  
        """