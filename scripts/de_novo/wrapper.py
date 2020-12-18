__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

SCRIPT_NAME = 'de_novo.py'

from snakemake.shell import shell
import os
import uuid
from tempfile import TemporaryDirectory, gettempdir

shell.executable('/bin/bash')

input_vcf = snakemake.input.get('vcf', None)
if input_vcf is None or not os.path.exists(input_vcf):
    raise Exception("Error: Please provide a valid input table.")

output_table = snakemake.output.get('table', None)
if output_table is None:
    raise Exception("Error: Please provide an output table path.")

fam = snakemake.params.get('fam', None)
if fam is None:
    raise Exception("Error: Please provide a pedigree file for de novo filter")

extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

memory = f"{snakemake.resources.get('mem_mb', 2000)}M"

script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), SCRIPT_NAME)
shell("""
    export JAVA_HOME="/home/brand/scratch/java/jdk8u265-b01"
    export PATH="/home/brand/scratch/java/jdk8u265-b01/bin:$PATH"
    export PYSPARK_SUBMIT_ARGS="--master local[{snakemake.threads}] --driver-memory {memory} pyspark-shell"
    python3 {script_path} --verbose debug {extra} --output {output_table} {input_vcf} {fam} {log}
""")