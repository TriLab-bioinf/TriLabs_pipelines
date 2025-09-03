# Call peaks with macs2
rule peak_calling:
    input:
        treat = "results/3-dedup/{sample}.dedup.bam",
    output:
        peaks = "results/4-peak_calling/{sample}_peaks.narrowPeak",
        control_bdg = "results/4-peak_calling/{sample}_control_lambda.bdg",
        chip_bdg = "results/4-peak_calling/{sample}_treat_pileup.bdg",
    params:
        gsize = config["reference"]["gsize"],
        name = "{sample}",
        outdir = "results/4-peak_calling"
    shell:
        """
        module load macs

        macs2 callpeak -t {input.treat} -f BAM -g {params.gsize} -n {params.name} -B -q 0.01 --outdir {params.outdir}

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
    shell:
        """
        module load ucsc
      
        cut -f1,2 {params.idx} > chrom.sizes
        bedGraphToBigWig {input.control} chrom.sizes {output.control_bw}
        bedGraphToBigWig {input.chip} chrom.sizes {output.chip_bw}
  
        """
