rule annotate:
    input:
        vcf='output/{sample}/{sample}.hard-filtered.genes.vcf.gz'
    output:
        vcf='output/{sample}/{sample}.hard-filtered.annotated.vcf.gz'
    params:
        extra=' '.join([
            '--sift', 's',
            '--polyphen', 's',
            # Gnomad exome
            '--af_gnomad',
            # Gnomad genome
            '--custom', ','.join([
                os.path.abspath(config['ref']['gnomad-genome']), 'gnomADg',
                'vcf', 'exact', '0', 'AF_AFR', 'AF_AMR', 'AF_ASJ', 'AF_EAS',
                'AF_FIN', 'AF_NFE', 'AF_OTH',
            ]),
            # CADD scores
            '--plugin', ','.join([
                'CADD',os.path.abspath(config['ref']['cadd']),
                os.path.abspath(config['ref']['cadd-indels'])
            ])
        ]),
        vep_cache=config['ref']['vep-cache'],
    log: "logs/annotate/{sample}.log"
    resources:
        runtime=60,
        mem_mb=8000,
    wrapper: "file:wrapper/vep"