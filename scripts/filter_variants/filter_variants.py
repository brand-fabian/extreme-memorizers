from __future__ import annotations
import pysam
import os
import pandas
import typing
import re
import argparse
import logging
import sys

##
# Logging
#
logger = logging.getLogger("filter-variants")
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

##
# Argument Parsing
#
parser = argparse.ArgumentParser(
    prog="FILTER-VARIANTS",
    description='Find variants in a given vcf file.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('vcf', metavar='VCF', type=str,
                    help='VCF input file for filtration.')
parser.add_argument('--annotation-field', default='CSQ', type=str,
                    help='Field containing VEP annotations')
parser.add_argument('-o', '--output', type=str, default=None,
                    help='Output file path', required=True)
parser.add_argument('--min-cadd', type=float, default=24,
                    help='Minimum cadd score for all variants')
parser.add_argument('--max-gnomad-af', type=float, default=0,
                    help='Maximum gnomad allele frequency')
parser.add_argument('--impact', nargs='+', default=[ 'HIGH', 'MODERATE' ],
                    help='Variant impact')

parser.add_argument('--verbose', '-v', help='Set verbosity',
                    choices=log_level.keys(), default='info')

args = parser.parse_args()
logger.setLevel(log_level[args.verbose])
logger.debug("Arguments: {}".format(vars(args)))

if not os.path.isfile(args.vcf):
    logger.error(f"Could not find vcf file {args.vcf}")
    sys.exit(1)

##
# Class and Function definitions
#
apply_fn = typing.Callable[[pysam.libcbcf.VariantRecord], pysam.libcbcf.VariantRecord]
filter_fn = typing.Callable[[pysam.libcbcf.VariantRecord], bool]
record_type = pysam.libcbcf.VariantRecord


record_attributes = {
    'alts': 'list',
    'chrom': 'fixed',
    'filter': 'list',
    'format': 'list',
    'id': 'fixed',
    'info': 'dict',
    'pos': 'fixed',
    'qual': 'fixed',
    'ref': 'fixed',
    'samples': 'dict',
}


class FilterIteratorException(Exception):
    pass


class RecordAttributeError(Exception):
    pass


class FilterVariants():
    class _TransformFn():
        def __init__(self, fn: typing.Union[apply_fn, filter_fn], is_filter: bool = True):
            self._fn = fn
            self.is_filter = is_filter
        
        def __call__(self, record: record_type):
            return self._fn(record)
        
    def __init__(self, vcf_file: str):
        self._vcf = pysam.VariantFile(vcf_file)
        self._transforms = []
        self._is_iterating = False
        self._filtered_records = 0
        self._annotation_header = None
        
    def apply(self, fn: apply_fn) -> FilterVariants:
        if self._is_iterating:
            raise FilterIteratorException("Changing transformations during iteration.")
        self._transforms.append(FilterVariants._TransformFn(fn, is_filter=False))
        return self
        
    def filter(self, fn: fiter_fn) -> FilterVariants:
        if self._is_iterating:
            raise FilterIteratorException("Changing transformations during iteration.")
        self._transforms.append(FilterVariants._TransformFn(fn, is_filter=True))
        return self
        
    def _apply_transforms(self, record: record_type) -> typing.Tuple[bool, record_type]:
        is_filtered = False
        for fn in self._transforms:
            if is_filtered:
                break
            if fn.is_filter:
                is_filtered = not fn(record)
            else:
                record = fn(record)
        return (is_filtered, record)
        
    @property
    def records(self):
        for record in self._vcf.fetch():
            self._is_iterating = True
            is_filtered, record = self._apply_transforms(record)
            if not is_filtered:
                yield record
            else:
                self._filtered_records += 1
        self._is_iterating = False
    
    @property
    def annotation_header(self) -> typing.List[str]:
        if self._annotation_header is None:
            self._annotation_header = list(
                filter(
                    lambda x: x.get('ID') == args.annotation_field,
                    self._vcf.header.records
                ))[0].get('Description').split(':')[1].strip(" \"").split('|')
        return self._annotation_header
    
    @classmethod
    def get_annotations(cls,
        record: record_type,
        header: typing.List[str],
        info_field=args.annotation_field
    ) -> typing.List[typing.Mapping[str, str]]:
        def _parse_val(val: str) -> typing.Union[str, float, int]:
            if re.match(r"^[+-]?\d+$", val):
                return int(val)
            elif re.match(r'^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$', val):
                return float(val)
            elif val == '':
                return None
            else:
                return val
        return [{
            header[i]: _parse_val(k) for i, k in enumerate(csq_field.split('|'))
        } for csq_field in record.info[info_field] ]
    
    def _sample_to_row(self,
        record: record_type,
        attributes: typing.List[str] = None,
    ) -> typing.Tuple[typing.List[typing.List[str]], typing.List[str]]:
        rows = []
        if attributes is None:
            s = list(record.samples)[0]
            s_attrs = set(list(record.samples[s]))
        else:
            attrs = set(filter(lambda s: s.startswith('samples'), attributes))
            s_attrs = set(map(lambda s: s.split('.')[1], attrs))
        header = None
        for s in list(record.samples):
            attrs_to_row = set(list(record.samples[s])) & s_attrs
            header = list(sorted(attrs_to_row))
            rows.append([
                s,
                *[ record.samples[s][a] for a in sorted(attrs_to_row) ]
            ])
        return rows, ['sample_id', *[ f'sample.{a}' for a in header ]]
            
    def _attrs_to_row(self,
        record: record_type,
        attributes: typing.List[str] = None,
        excl_attr: typing.List[str] = [ 'samples' ],
        anno_field: str = args.annotation_field,
        split_anno: bool = True,
    ) -> typing.Tuple[typing.List[typing.List[str]], typing.List[str]]:
        
        def _get_single_row(record, attrs, anno_field=anno_field) -> typing.List[str]:
            row = []
            for a in attrs:
                a_key = a.split('.')[0]
                if record_attributes[a_key] == 'fixed':
                    row.append(getattr(record, a_key))
                elif record_attributes[a_key] == 'list':
                    row.append(','.join(list(getattr(record, a_key))))
                elif record_attributes[a_key] == 'dict':
                    attr_value = getattr(record, a_key)
                    if len(a.split('.')) == 1:
                        # Get all
                        for key in list(attr_value):
                            if a_key == 'info' and key == anno_field:
                                continue
                            else:
                                row.append(attr_value[key])
                    elif a.split('.')[1] in list(attr_value):
                        row.append(attr_value[a.split('.')[1]])
                    else:
                        raise RecordAttributeError(f"{a} could not be found on record.")
            return row

        if attributes is not None:
            attrs = set(filter(lambda s: s.split('.')[0] in record_attributes.keys())) - set(excl_attr)
        else:
            attrs = set(record_attributes.keys()) - set(excl_attr)

        has_anno = any(a.startswith(f'info.{anno_field}') for a in attrs)
        if has_anno or attributes is None:
            anno = FilterVariants.get_annotations(record, self.annotation_header)
        else:
            anno = None
            
        if anno is None:
            return [ _get_single_row(record, attrs) ]
        else:
            if len(set(list(record.info)) - set([ anno_field ])) == 0:
                attrs -= set([ 'info' ])
            return [
                [ *a.values(), *_get_single_row(record, attrs) ] for a in anno 
            ], [ *[ f"info.{a}" for a in anno[0].keys() ] , *attrs ]

    
    def _record_to_row(self,
        record: record_type,
        attributes: typing.List[str] = None,
        info_field: str = args.annotation_field,
        split_info: bool = True
    ) -> typing.Tuple[typing.List[typing.List[str]], typing.List[str]]:
        rows = []
        sample_rows = self._sample_to_row(record, attributes)
        row_attrs = self._attrs_to_row(record, attributes, anno_field=info_field, split_anno=split_info)
        var_id = f"{record.contig}_{record.start}_{record.ref}_{record.alts[0]}"
        header = [ 'var_id', *row_attrs[1], *sample_rows[1] ]
        for s in sample_rows[0]:
            for a in row_attrs[0]:
                rows.append([
                    var_id,
                    *a,
                    *s,
                ])
        return rows, header 
    
    def to_pandas(self, attributes: typing.List[str] = None) -> pandas.DataFrame:
        rows = []
        header = []
        for record in self.records:
            data, header = self._record_to_row(record, attributes)
            rows = [ *rows, *data ]
        return pandas.DataFrame(rows, columns=header)


##
# Main script
#

logger.info("Filter variants...")
filter_obj = FilterVariants(args.vcf)


def filter_cadd(
    record: record_type,
    header: typing.List[str] = filter_obj.annotation_header,
    min_cadd: int = args.min_cadd
) -> bool:
    anno = FilterVariants.get_annotations(record, header)
    return any('CADD_PHRED' in a.keys() and a['CADD_PHRED'] is not None and a['CADD_PHRED'] >= min_cadd for a in anno)

def filter_gnomad(
    record: record_type,
    header: typing.List[str] = filter_obj.annotation_header,
    max_gnomad_af: float = args.max_gnomad_af
) -> bool:
    anno = FilterVariants.get_annotations(record, header)
    return any('gnomAD_AF' in a and (a['gnomAD_AF'] is None or a['gnomAD_AF'] <= max_gnomad_af) for a in anno)

def filter_impact(
    record: record_type,
    header: typing.List[str] = filter_obj.annotation_header,
    impact: typing.List[str] = args.impact
) -> bool:
    anno = FilterVariants.get_annotations(record, header)
    return any('IMPACT' in a and a['IMPACT'] is not None and a['IMPACT'] in impact for a in anno)

filter_obj = filter_obj.filter(filter_impact).filter(filter_gnomad).filter(filter_cadd)
df = filter_obj.to_pandas()

logger.info(f"Filtered {filter_obj._filtered_records} records.")
if len(df) > 0 and 'info.IMPACT' in df.columns:
    df = df[df['info.IMPACT'].apply(lambda s: s in args.impact)]

logger.info(f"Writing {len(df)} records to {args.output}")
df.to_csv(args.output, sep="\t", index=False, header=True)

logger.info("Done.")