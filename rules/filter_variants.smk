rule filter_variants:
    input:
        vcf='output/{sample}/{sample}.hard-filtered.annotated.vcf.gz'
    output:
        table='output/{sample}/{sample}.variants.tsv'
    params:
        extra=' '.join([
            '--lof-metrics', config['ref']['lof-metrics']
        ])
    log: 'log/filter_variants/{sample}.log'
    resources:
        runtime=60,
        mem_mb=8000,
    wrapper: 'file:scripts/filter_variants'