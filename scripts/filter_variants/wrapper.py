__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

SCRIPT_NAME = 'filter_variants.py'

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


extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), SCRIPT_NAME)
shell("""
    python3 {script_path} --verbose debug --output {output_table} {input_vcf} {log}
""")