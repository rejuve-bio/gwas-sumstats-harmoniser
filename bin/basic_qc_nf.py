#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import logging

import pandas as pd

from common_constants import *

logger = logging.getLogger('basic_qc')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(message)s')

hm_header_transformations = {
    'hm_varid': HM_VAR_ID,
    'hm_OR': HM_OR_DSET,
    'hm_OR_lowerCI': HM_RANGE_L_DSET,
    'hm_OR_upperCI': HM_RANGE_U_DSET,
    'hm_beta': HM_BETA_DSET,
    'hm_effect_allele': HM_EFFECT_DSET,
    'hm_other_allele': HM_OTHER_DSET,
    'hm_eaf': HM_FREQ_DSET,
    'hm_code': HM_CODE
}

REQUIRED_HEADERS = [RSID, PVAL_DSET, CHR_DSET, BP_DSET]
BLANK_SET = {'', ' ', '-', '.', 'na', None, 'none', 'nan', 'nil'}

# hm codes to drop (kept as ints to preserve legacy type-mismatch behaviour —
# string hm_code values in the file will NOT match these ints, so this filter
# is effectively a no-op on the data, matching the original csv-reader version)
HM_CODE_FILTER = {9, 14, 15, 16, 17, 18}


def check_for_required_headers(header):
    return list(set(REQUIRED_HEADERS) - set(header))


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-f', help='Input harmonised file', required=True)
    argparser.add_argument('-o', help='Output file', required=True)
    argparser.add_argument('-db', help='Synonyms database (unused, kept for compatibility)', default=None)
    argparser.add_argument('--print_only', help='Log removed rows without writing output', action='store_true')
    argparser.add_argument('--log', help='Log file path')
    args = argparser.parse_args()

    file_handler = logging.FileHandler(args.log, mode='a')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Read entire file with pandas (much faster than csv.reader row-by-row)
    df = pd.read_csv(args.f, sep='\t', dtype=str)

    # Rename legacy column headers
    df.columns = [hm_header_transformations.get(c, c) for c in df.columns]

    missing_headers = check_for_required_headers(list(df.columns))
    if missing_headers:
        sys.exit("Required headers are missing!:{}".format(missing_headers))

    # Replace blank-like values with 'NA' across all columns (vectorized)
    blank_lower = {v.lower() for v in BLANK_SET if v is not None and isinstance(v, str)}
    blank_lower.add('')
    df = df.apply(
        lambda col: col.where(~col.str.lower().isin(blank_lower), other='NA')
    )

    # Map chromosome strings X/Y/MT → 23/24/25
    chr_lower = df[CHR_DSET].str.lower()
    chr_map = {'x': '23', 'y': '24', 'mt': '25'}
    df[CHR_DSET] = chr_lower.map(lambda x: chr_map.get(x, x) if isinstance(x, str) else x)

    # Determine which rows to remove
    # Note: HM_CODE_FILTER contains ints; df[HM_CODE] contains strings.
    # The isin comparison intentionally does not match (legacy behaviour preserved).
    mask_unharmonisable = df[HM_CODE].isin(HM_CODE_FILTER)

    na_values = {'NA', 'na', 'none', 'nan', 'nil', '', ' ', '-', '.'}
    mask_blank_required = df[REQUIRED_HEADERS].apply(
        lambda col: col.str.lower().isin(na_values)
    ).any(axis=1)

    mask_bad_chr  = pd.to_numeric(df[CHR_DSET],  errors='coerce').isna()
    mask_bad_bp   = pd.to_numeric(df[BP_DSET],   errors='coerce').isna()
    mask_bad_pval = pd.to_numeric(df[PVAL_DSET], errors='coerce').isna()

    mask_remove = mask_unharmonisable | mask_blank_required | mask_bad_chr | mask_bad_bp | mask_bad_pval

    # Log removed rows
    for idx in df.index[mask_remove]:
        try:
            hm_code = int(df.loc[idx, HM_CODE])
        except (ValueError, TypeError):
            hm_code = 19
        if hm_code not in HM_CODE_FILTER:
            hm_code = 19
        logger.info(f'Removing record number {idx + 1}, with hm_code {hm_code}')

    if not args.print_only:
        df[~mask_remove].to_csv(args.o, sep='\t', index=False)


if __name__ == "__main__":
    main()
