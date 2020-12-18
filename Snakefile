include: "rules/meta.smk"

##
# Target files
#
rule all:
    input:
        expand("output/{sample}/{sample}.hard-filtered.genes.vcf.gz", sample=samples.keys()),
        expand("output/{sample}/{sample}.hard-filtered.annotated.vcf.gz", sample=samples.keys()),
        expand("output/{sample}/{sample}.variants.tsv", sample=samples.keys()),
        expand("output/{sample}/{sample}.de_novo.csv", sample=samples.keys()),

##
# Snakemake rules
#
include: "rules/extract_genes.smk"
include: "rules/annotate.smk"
include: "rules/filter_variants.smk"
include: "rules/de_novo.smk"