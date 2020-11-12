rule extract_genes:
    input:
        vcf=get_vcf
    output:
        vcf="output/{sample}/{sample}.hard-filtered.genes.vcf.gz"
    params:
        extra=" ".join([
            "--remove-filtered-all",
            "--bed",
            config['ref']['covered-regions'],
        ])
    log: "logs/extract_genes/{sample}.log"
    resources:
        runtime=120,
        mem_mb=16000,
    wrapper: "file:wrapper/vcftools"