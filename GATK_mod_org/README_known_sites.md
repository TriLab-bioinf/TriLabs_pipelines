# Using known SNP/INDEL sites

## Which Known-Sites File Should You Use for _Drosophila melanogaster_?

For _Drosophila melanogaster_, you can run GATK BaseRecalibrator, but the key challenge is having a trusted set of known variants so that true biological variants are not mistaken for sequencing errors.

Below are the best available and commonly used known-sites resources for D. melanogaster.

### Recommended Known-Sites Databases for _Drosophila melanogaster_

Unlike human data, where dbSNP and other curated resources exist, fly genomics relies mainly on population projects.

1. **Drosophila Genetic Reference Panel (DGRP)**

    - Best and most widely used source for known SNPs/indels

    - High-quality variant calls from ∼200 inbred lines

    - Genome version: usually dm3 (BDGP release 5), but many liftover resources exist to dm6

    Download:

    Website: https://dgrp2.gnets.ncsu.edu/data.html

    Direct VCFs also available from NCBI SRA or from prior publications’ supplementary files.

    These variants are commonly used as the “known sites” input for BQSR.

2. **Drosophila Population Genomics Project (DPGP)**

    - Covers African and cosmopolitan populations.

    - Provides SNP and indel datasets for several collections.

    - Good for representing natural variation.

    Links:

    DPGP2 & DPGP3 projects (VCFs available)

3. **FlyBase Variation Data**

    - Contains curated SNP/indel annotations.

    - Derived from various sequencing projects, including modENCODE.

    - Can be used as a general known-variants VCF.

    Download:

    https://flybase.org

    Under: Genomes → Variations → Polymorphisms (SNPs/Indels)

    Note: Data may not be as population-specific as DGRP.

4. **ModENCODE Project Variant Calls**

    - Large comparative genomics effort.

    - Has SNP/indel datasets for specific inbred strains.

    - Useful if you are working with corresponding strains.

#### Which One Should You Use?

Best practice (similar to human GATK pipelines):

- If working with lab strains related to DGRP (very common):

    Use DGRP SNPs + DGRP indels

- If working with wild-caught or diverse natural populations:

    Combine both: DPGP SNP/indels + FlyBase polymorphisms

- If working with reference strains such as ISO1:

    Use FlyBase variations or curated ISO1 variant sets.

⚠ **Important Notes About BQSR in Non-Human Genomes**

1. A perfect known-set is not required

    GATK recommends:

    A large set of reasonably accurate variants is enough.

    False positives are tolerated as long as they are not systematic.

    This means using DGRP or DPGP variant sets works well even if not exact to your strain.

2. The more variation your sample has relative to reference, the more crucial the known-sites file becomes.
3. If no known-sites are available at all, you can:

    Generate initial variant calls

    Filter for high-confidence SNPs

    Use those as known sites

    Re-run BQSR

    This is known as bootstrapping (similar to tumor-only pipelines).

    _Example: How to run BaseRecalibrator with DGRP known sites_

```bash
gatk BaseRecalibrator \
    -I sample.bam \
    -R dmel_reference.fasta \
    --known-sites dgrp2_snps.vcf.gz \
    --known-sites dgrp2_indels.vcf.gz \
    -O recal_data.table
```

## Which Known-Sites File Should You Use for _Caenorhabditis elegans_?

For Caenorhabditis elegans (C. elegans), GATK BaseRecalibrator can be used, but as with other non-human organisms, the challenge is obtaining a good known-variants VCF so that true polymorphisms are not counted as sequencing errors.

Fortunately, several high-quality variation resources exist for C. elegans.

Below are the best and most commonly used known-sites datasets for BQSR in C. elegans.

### Recommended Known-Sites Databases for _C. elegans_

1. CeNDR – Caenorhabditis elegans Natural Diversity Resource

    ➡ Best and most widely used choice for known SNP/indel sites
    CeNDR provides high-quality variant calls from hundreds of C. elegans wild isolates, aligned to the WBcel235 (same as ce11) reference.

    What CeNDR offers:

    SNPs and indels in VCF format

    Hard-filtered, high-confidence calls

    Population-wide allele frequencies

    Consistent reference version (WS245 / WBcel235)

    Download (VCFs):

    https://www.elegansvariation.org/data/release

    Typical files:

    WI.20170531.hard-filter.vcf.gz
    WI.20200815.hard-filter.vcf.gz

    These are ideal for GATK BQSR.

2. WormBase Variation Data

    Alternative or additional known-sites set.

    WormBase provides:

    Curated variants

    Known polymorphisms

    Strain-specific mutations (e.g., Hawaiian CB4856 strain)

    Download:

    https://wormbase.org
    → Tools → Variation

    VCFs can be downloaded for specific collections.

    These datasets are good if:

    You work with the canonical N2 strain

    You want curated variants instead of population variation

3. CB4856 / Hawaiian Strain Variants

    If your sample is closely related to CB4856 (one of the most divergent strains), you may want to include known variants specific to that lineage.

    Sources:

    CeNDR (strain-specific VCFs)

    Various population genetics publications

    Useful when sequencing crosses or introgression lines involving CB4856.

4. Mutation Accumulation (MA) Lines Variant Sets

    For experiments involving:

    Mutator strains

    Mutation accumulation studies

    You can use known variants from specific MA line studies.

    Common studies:

    Denver et al. MA lines

    Katju & Bergthorsson MA datasets

    But these are not recommended for general BQSR—only for specialized experiments.

## Which Known-Sites File Should You Use?

1. **Best all-purpose choice (recommended):**

    CeNDR hard-filtered SNP + indel VCF

2. **For wild isolates or general natural diversity studies:**

    CeNDR + WormBase variation

3. **For laboratory N2 strain whole-genome sequencing:**

    WormBase curated variants (N2-specific)

4. **For CB4856 or hybrid/mapping populations:**

    Add CB4856-specific variant VCFs (from CeNDR)

If you also use indels:

--known-sites WI.20200815.hard-filter.indels.vcf.gz

## Which Known-Sites File Should You Use for _Escherichia coli_?

Using GATK BaseRecalibrator for Escherichia coli (_E. coli_) is possible, but unlike human, fly, or worm genomes, there is no widely used, standardized “known variants” VCF for _E. coli_. This is because _E. coli_ is highly diverse, strains differ substantially, and there is no single population reference panel.

However, you can still run BaseRecalibrator using appropriate known‐sites resources or by bootstrapping.

Below are the recommended approaches.

1. **If you are sequencing a well-studied lab strain (K-12 MG1655 or BW25113)**

    There are published variant sets that are commonly used, especially for the classical laboratory strain MG1655.

    Option A — Curated MG1655 mutation datasets

    Several studies and databases list:

    long-term evolution experiment mutations

    verified MG1655 polymorphisms

    other curated SNP/indel sets

    Examples:

    E. coli K-12 reference corrections from Latif et al., 2014

    Breseq reference corrections for MG1655

    KEGG/RefSeq MG1655 annotation discrepancies

    These can be merged into a known-sites VCF.

    If you want, I can generate a clean VCF for MG1655 from these curated sources.

2. **If you are sequencing a strain closely related to MG1655 or BW25113**

    Use known variants from:

    EcoCyc / EcoGene corrections

    Breseq reference corrections
    (Breseq provides reference "fix files" for MG1655 and BW25113)

    Example:
    https://barricklab.org/breseq

    These are typically small but high-confidence sets of corrections to the reference.

    Not population-wide, but acceptable for BQSR.

3. **If you are sequencing a non-classical strain (pathogenic or environmental)**

    There is no population-wide variant set comparable to dbSNP.

    Your best option is bootstrap BQSR (recommended for bacteria)

    This is the standard approach for microbial genomes:

    - Steps 1:

        Align reads

        Call variants (e.g., using HaplotypeCaller or FreeBayes)

        Filter to very high-confidence SNPs
        (QUAL > 100, DP > 20, no strand bias)

    - Step 2

        Use the filtered SNPs as --known-sites for BaseRecalibrator

        Re-run variant calling

    - Step 3 (optional)

        repeat once more

    This produces very good recalibration without needing external databases. Bootstrapping is simple and robust for microbial genomes because:

    - Genome is small

    - Coverage is usually high

    - True variation is consistent within a strain

4. **If you have multiple isolates of the same strain/species**

    You can construct a pan-genome SNP set:

    - Call variants across all isolates

    - Keep only shared, high-support SNPs

    - Use this as a known-sites VCF for BQSR

    This works especially well for:

    - Evolution experiments

    - Clinical sampling

    - Environmental metapopulation sampling

❌ **What you should NOT use**

Generic _E. coli_ SNP sets from NCBI/ENA:

- These are not compatible with your specific strain's reference genome.

Pan-genome VCFs from broad _E. coli_ collections (e.g., 50,000 isolate catalogs).

- Too many sites, wrong coordinates, not strain-specific.

These will harm BQSR more than help.

**In summary:**

a- Sequencing MG1655 or BW25113:

Use curated MG1655/BW25113 correction SNP sets (Breseq, EcoCyc, published corrections)

b- Sequencing strains similar to K-12

Same as above; optional bootstrap

c- Sequencing any other _E. coli_ strain

Bootstrap known-sites from your own data

d- Many isolates from same lineage

Build a high-confidence shared SNP VCF

_Example of bootstrap approach:_

Round 1: Call variants

```bash
gatk HaplotypeCaller -I sample.bam -R ecoli.fasta -O raw.vcf
```

Filter high-confidence SNPs

```bash
gatk VariantFiltration \
  -V raw.vcf \
  --filter-expression "QD < 10.0 || FS > 10.0 || MQ < 40.0" \
  --filter-name "FAIL" \
  -O filtered.vcf
```

Keep only PASS sites:

```bash
gatk SelectVariants -V filtered.vcf --exclude-filtered -O known_sites.vcf
```

Run BaseRecalibrator

```bash
gatk BaseRecalibrator \
  -I sample.bam \
  -R ecoli.fasta \
  --known-sites known_sites.vcf \
  -O recal_data.table
```
