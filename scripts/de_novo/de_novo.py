import hail as hl
import logging
import argparse
import os
import sys

##
# Logging
logger = logging.getLogger("de-novo")
logger.setLevel(logging.ERROR)
channel = logging.StreamHandler()
formatter = logging.Formatter("[%(asctime)s - %(name)s - %(levelname)7s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
channel.setFormatter(formatter)
logger.addHandler(channel)

log_level = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
}

output_types = {
    'PANDAS': lambda ht, path: ht.to_pandas().to_csv(path),
    'HAIL': lambda ht, path: ht.write(path, overwrite=True),
}

parser = argparse.ArgumentParser(
    prog="DE-NOVO",
    description='Find de-novo variants in a given vcf file.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('vcf', metavar='VCF', type=str,
                    help='Trio VCF input file for filtration.')
parser.add_argument('fam', metavar='FAM', type=str,
                    help='.fam pedigree file.')
parser.add_argument('-o', '--output', type=str, default=None,
                    help='Output file path', required=True)
parser.add_argument('-O', '--output-type', choices=output_types.keys(),
                    help='Output data type', default='PANDAS')
parser.add_argument('--min-cadd', type=float, default=24,
                    help='Minimum cadd score for all variants')
parser.add_argument('--max-gnomad-af', type=float, default=0,
                    help='Maximum gnomad allele frequency')
parser.add_argument('--impact', nargs='+', default=[ 'HIGH', 'MODERATE' ],
                    help='Variant impact')
parser.add_argument('--lof-metrics', type=str, default=None,
                    help='Gnomad lof metrics table')
parser.add_argument('--min-pli', type=float, default=0.9,
                    help='Min. pli prob. for lof genes')
parser.add_argument('--max-oe-lof', type=float, default=0.35,
                    help='Max. observed expected ratio for lof genes')
parser.add_argument('--min-mis-z', type=float, default=3.09,
                    help='Min. missense mutation z score')
parser.add_argument('--min-p-de-novo', type=float, default=0.80,
                    help='Min. probability of de novo event')
parser.add_argument('--min-aaf', type=float, default=0.30,
                    help='Min. alternate allele frequency to call het')
parser.add_argument('--min-parent-dp', type=int, default=10,
                    help='Min. sequencing depth in parents')
parser.add_argument('--min-dp', type=int, default=15,
                    help='Min. sequencing depth in index')
parser.add_argument('--max-parent-alleles', type=int, default=1,
                    help='Max. reads in parents supporting alt allele')
parser.add_argument('--missing-true', default=True, type=bool,
                    help='Treat missing values as true in filter settings')
parser.add_argument('-t', '--tmp-dir', help='Temporary directory',
                    type=str, default=os.getenv('TMPDIR', '/tmp'))
parser.add_argument('--dnm-checkpoint', help='Write dnm checkpoint data',
                    type=str, default=None)
parser.add_argument('--vep-checkpoint', help='Write vep checkpoint data',
                    type=str, default=None)
parser.add_argument('--vep-config', default='./vep_config.json',
                    type=str, help='Path to hail vep configuration')
parser.add_argument('--verbose', '-v', help='Set verbosity',
                    choices=log_level.keys(), default='info')

args = parser.parse_args()
logger.setLevel(log_level[args.verbose])
logger.debug("Arguments: {}".format(vars(args)))

if not os.path.isfile(args.vcf):
    logger.error(f"Could not find vcf file {args.vcf}.")
    sys.exit(1)

if not os.path.isfile(args.fam):
    logger.error(f"Could not find pedigree file {args.fam}.")
    sys.exit(1)

if args.lof_metrics is not None and not os.path.isfile(args.lof_metrics):
    logger.error(f"Could not read gnomad lof metrics file {args.lof_metrics}.")
    sys.exit(1)
else:
    filter_constraint = True

if not os.path.isfile(args.vep_config):
    logger.error(f"Could not find vep config file {args.vep_config}")
    sys.exit(1)


# Prepare output path
if not os.path.exists(os.path.abspath(os.path.dirname(args.output))):
    os.makedirs(os.path.abspath(os.path.dirname(args.output)))

# Set hail temporary path
hl.init(
    idempotent=True,
    tmp_dir=args.tmp_dir,
    log=os.path.join(args.tmp_dir, 'hail.log'),
)


##
# Main script
#
logger.info(f"Reading pedigree file {args.fam}")
pedigree = hl.Pedigree.read(args.fam)

logger.info(f"Importing vcf file {args.vcf}")
data = hl.import_vcf(args.vcf, call_fields=['GT'],
                     skip_invalid_loci=True, force_bgz=True)
data = hl.split_multi_hts(data)
data = data.annotate_rows(
    AC=data.info.AC[data.a_index - 1],
    iAF=data.info.AF[data.a_index - 1]
)
data = hl.variant_qc(data)

logger.info("Applying de novo filter...")
de_novo_scores = hl.de_novo(data, pedigree, pop_frequency_prior=data.variant_qc.AF[-1])
de_novo_mt = de_novo_scores.to_matrix_table(row_key=['locus', 'alleles'], col_key=['id'])
de_novo_data = data.annotate_entries(p_de_novo=de_novo_mt[(data.locus, data.alleles), data.s].p_de_novo)

logger.info("Annotating trio data...")
trio_mt = hl.trio_matrix(de_novo_data, pedigree, complete_trios=True)
de_novo_data = de_novo_data.annotate_entries(
    mother=trio_mt[(de_novo_data.locus, de_novo_data.alleles), de_novo_data.s].mother_entry,
    father=trio_mt[(de_novo_data.locus, de_novo_data.alleles), de_novo_data.s].father_entry,
)
de_novo_data = de_novo_data.filter_entries(hl.is_defined(de_novo_data.GT)
    & hl.is_defined(de_novo_data.PL)
    & de_novo_data.GT.is_non_ref()
    & (de_novo_data.p_de_novo > args.min_p_de_novo)
)

r_de_novo_mt = de_novo_data.select_cols()
r_de_novo_mt = r_de_novo_mt.select_rows('AC', 'iAF')
r_de_novo_mt = r_de_novo_mt.select_entries('AD', 'AF', 'DP', 'DN', 'DQ', 'GT', 'p_de_novo',
    mother=hl.Struct(**{
        'AD': r_de_novo_mt.mother.AD,
        'AF': r_de_novo_mt.mother.AF,
        'DP': r_de_novo_mt.mother.DP,
        'GT': r_de_novo_mt.mother.GT
    }),
    father=hl.Struct(**{
        'AD': r_de_novo_mt.father.AD,
        'AF': r_de_novo_mt.father.AF,
        'DP': r_de_novo_mt.father.DP,
        'GT': r_de_novo_mt.father.GT        
    })
)

logger.info("Refining de-novo calls...")
if args.dnm_checkpoint is None or not os.path.exists(args.dnm_checkpoint):
    dnm = r_de_novo_mt.filter_entries(
        (r_de_novo_mt.p_de_novo > args.min_p_de_novo)
        & (r_de_novo_mt.DP >= args.min_dp)
        & (r_de_novo_mt.mother.DP >= args.min_parent_dp)
        & (r_de_novo_mt.father.DP >= args.min_parent_dp)
        & (r_de_novo_mt.AF[0] > args.min_aaf)
        & (r_de_novo_mt.AF[0] < 1 - args.min_aaf)
    )
    dnm_ht = dnm.key_cols_by().entries()
    dnm_ht = dnm_ht.filter(
        (dnm_ht.alleles[0].length()==1)
        & (dnm_ht.alleles[1].length()==1)
    )
    dnm_ht = dnm_ht.filter(
        (dnm_ht.father.AD[1] <= args.max_parent_alleles)
        & (dnm_ht.mother.AD[1] <= args.max_parent_alleles)
    )
    if args.dnm_checkpoint is not None:
        dnm_ht = dnm_ht.checkpoint(args.dnm_checkpoint)
else:
    dnm_ht = hl.read_table(args.dnm_checkpoint)

logger.info("Annotating data using vep...")
if args.vep_checkpoint is None or not os.path.exists(args.vep_checkpoint):
    vep_dnm_ht = hl.vep(dnm_ht, config=args.vep_config, name="vep", csq=False)
    if args.vep_checkpoint is not None:
        vep_dnm_ht = vep_dnm_ht.checkpoint(args.vep_checkpoint, overwrite=True)
else:
    vep_dnm_ht = hl.read_table(args.vep_checkpoint)

def max_if_defined(expr):
    return hl.nanmax(hl.if_else(hl.is_defined(expr),
        expr,
        [ .0 ],
    ))

vep_dnm_ht = vep_dnm_ht.annotate(
    gene_symbol=hl.set(vep_dnm_ht.vep.transcript_consequences.gene_symbol),
    gnomad=hl.Struct(**{
        'POPMAX': hl.nanmax([
            max_if_defined(vep_dnm_ht.vep.custom_annotations.gnomADg.fields.AF_AFR),
            max_if_defined(vep_dnm_ht.vep.custom_annotations.gnomADg.fields.AF_AMR),
            max_if_defined(vep_dnm_ht.vep.custom_annotations.gnomADg.fields.AF_ASJ),
            max_if_defined(vep_dnm_ht.vep.custom_annotations.gnomADg.fields.AF_EAS),
            max_if_defined(vep_dnm_ht.vep.custom_annotations.gnomADg.fields.AF_FIN),
            max_if_defined(vep_dnm_ht.vep.custom_annotations.gnomADg.fields.AF_NFE),
            max_if_defined(vep_dnm_ht.vep.custom_annotations.gnomADg.fields.AF_OTH),
        ]),
        'FILTER': vep_dnm_ht.vep.custom_annotations.gnomADg.fields.FILTER,
    }),
    cadd=hl.Struct(**{
        'is_defined': (
            hl.is_defined(vep_dnm_ht.vep.intergenic_consequences.cadd_phred)
            | hl.is_defined(vep_dnm_ht.vep.motif_feature_consequences.cadd_phred)
            | hl.is_defined(vep_dnm_ht.vep.regulatory_feature_consequences.cadd_phred)
            | hl.is_defined(vep_dnm_ht.vep.transcript_consequences.cadd_phred)
        ),
        'phred': hl.nanmax([
            max_if_defined(vep_dnm_ht.vep.intergenic_consequences.cadd_phred),
            max_if_defined(vep_dnm_ht.vep.motif_feature_consequences.cadd_phred),
            max_if_defined(vep_dnm_ht.vep.regulatory_feature_consequences.cadd_phred),
            max_if_defined(vep_dnm_ht.vep.transcript_consequences.cadd_phred),
        ]),
        'raw': hl.nanmax([
            max_if_defined(vep_dnm_ht.vep.intergenic_consequences.cadd_raw),
            max_if_defined(vep_dnm_ht.vep.motif_feature_consequences.cadd_raw),
            max_if_defined(vep_dnm_ht.vep.regulatory_feature_consequences.cadd_raw),
            max_if_defined(vep_dnm_ht.vep.transcript_consequences.cadd_raw),
        ]),
    }),
    impact=hl.set(vep_dnm_ht.vep.intergenic_consequences.impact)\
        .union(hl.set(vep_dnm_ht.vep.motif_feature_consequences.impact))\
        .union(hl.set(vep_dnm_ht.vep.regulatory_feature_consequences.impact))\
        .union(hl.set(vep_dnm_ht.vep.transcript_consequences.impact)),
    consequences=vep_dnm_ht.vep.most_severe_consequence
)
vep_dnm_ht = vep_dnm_ht.key_by().select(
    'locus', 'alleles', 'AC', 'iAF', 's', 'AD', 'AF', 'DP', 'GT', 'p_de_novo',
    'mother', 'father', 'gene_symbol', 'gnomad', 'cadd', 'impact', 'consequences',
)

logger.info("Filtering data with cadd and gnomad criteria...")
results_ht = vep_dnm_ht.filter(
    hl.if_else(vep_dnm_ht.cadd.is_defined,
        vep_dnm_ht.cadd.phred > args.min_cadd,
        args.missing_true
    )
    & (vep_dnm_ht.gnomad.POPMAX <= args.max_gnomad_af)
)

logger.info("Writing results...")
output_types[args.output_type](results_ht, args.output)