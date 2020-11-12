from snakemake.utils import validate, min_version
from snakemake.logging import logger
from pathlib import Path
from collections import namedtuple
import itertools
import json
import os

##
# Snakemake setup
#
report: "../report/extreme-memorizers.rst"

##
# Configuration files
#
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

##
# Setup scratch directory handling
#
work_dir = os.getenv('SCRATCH_DIR', os.getcwd())
logger.info(f"Working directory: {work_dir}")

##
# Helper functions
#
def find_files(base_path, file_endings=[], file_patterns=[]):
    base = Path(base_path)
    found_files = []
    for child in base.iterdir():
        if child.is_file():
            glob_patterns = [*file_patterns, *["*.{}".format(ending) for ending in file_endings]]
            if any(child.match(g) for g in glob_patterns):
                found_files.append(child)
                break
        if child.is_dir():
            found_files = [
                *found_files,
                *find_files(str(child), file_endings=file_endings,
                            file_patterns=file_patterns),
            ]
    return found_files

##
# Input files and wildcards setup
#
endings = config.get('file-endings', [])
patterns = config.get('search-patterns', [])
files = list(*itertools.chain(
    find_files(p, file_endings=endings,
               file_patterns=patterns) for p in config.get('input-dirs', [])
))
logger.info(f"Found {len(files)} files.")

samples = {}
separator = config.get('path-separator', '.')
for f in files:
    sample_name = f.name.split(separator)[0]
    #
    # Testing: Subset samples to one
    #
    if sample_name != '102-08068':
        continue
    #
    # End testing
    #
    if not sample_name in samples.keys():
        samples[sample_name] = set()
    samples[sample_name].add(str(f))
logger.info(f"Found {len(samples.keys())} samples.")

wildcard_constraints:
    sample="|".join(samples.keys()),
    project=config.get('project', 'project')

##
# Check required library files
#
def assert_path(path, is_file=True, name='path'):
    if path is None:
        logger.error(f"'{name}' is not set.")
        raise Exception(f"'{name}' is not set.")
    if is_file:
        if not os.path.isfile(path):
            logger.error(f"'{name}': file {path} does not exist.")
            raise Exception(f"'{name}': file {path} does not exist.")
    else:
        if not os.path.exists(path):
            logger.error(f"'{name}': path {path} does not exist.")
            raise Exception(f"'{name}': path {path} does not exist.")
ref = config.get('ref', None)
if ref is None:
    logger.error("'ref' configuration is missing. Please provide necessary files.")
    raise Exception("'ref' is missing.")
assert_path(ref.get('covered-regions', None), name='covered-regions')
assert_path(ref.get('vep-cache', None), is_file=False, name='vep-cache')
assert_path(ref.get('reference', None), is_file=True, name='reference')
assert_path(ref.get('genes', None), name='genes')
assert_path(ref.get('gnomad-exome', None), name='gnomad-exome')
assert_path(ref.get('gnomad-genome', None), name='gnomad-genome')
assert_path(ref.get('cadd', None), name='cadd')

##
# Input functions
#
def get_vcf(wildcards):
    vcfs = list(f for f in samples[wildcards.sample] if 'vcf' in f)
    if len(vcfs) > 1:
        logger.warning(f"WARNING: Found {len(vcfs)} vcfs for sample {wildcards.sample}")
    return vcfs[0]