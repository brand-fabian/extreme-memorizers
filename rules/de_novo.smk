rule de_novo:
    input:
        vcf=get_vcf
    output:
        table='output/{sample}/{sample}.de_novo.csv'
    params:
        fam=config['pedigree'],
        extra=' '.join([
            '--vep-config', config['vep-config'],
            '--dnm-checkpoint', 'output/{sample}/{sample}.dnm.ht',
            '--vep-checkpoint', 'output/{sample}/{sample}.vep.ht',
        ])
    log: 'logs/de_novo/{sample}.log'
    threads: 20
    resources:
        runtime=240,
        mem_mb=48000,
    wrapper: 'file:scripts/de_novo'
