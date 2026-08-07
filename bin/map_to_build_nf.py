#!/usr/bin/env python
# -*- coding: utf-8 -*-

# merge on hm_variant_id with vcf of desired build
# if variant_id is rsid and ID != variant_id or not synonym --> drop and create discrep df


from functools import lru_cache
import pandas as pd
import pyarrow.dataset as ds

import liftover as lft
from common_constants import *
import os
import glob
import argparse
from ast import literal_eval

# Allow very large fields in input file-------------
import sys
import csv

maxInt = sys.maxsize

while True:
    # decrease the maxInt value by factor 10
    # as long as the OverflowError occurs.

    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

# map_to_build----------------------------------------------------
# process the chr: if it is 23,24,25, convert to X,Y,MT
def normalize_chrom(c):
    return {"23": "X", "24": "Y", "25": "MT"}.get(str(c).upper(), str(c).upper())

def merge_ss_vcf(ss, vcf, from_build, to_build, chroms, coordinate, threads=1, memory="4GB"):

    """
    Merge GWAS summary stats with reference VCFs by RSID, liftover unmapped variants,
    and write per-chromosome outputs.

    Parameters:
    - ss (str): Path to summary statistics file
    - vcf (str): Glob pattern or exact path for reference VCF Parquet file(s)
    - from_build (str): Genome build of the summary stats
    - to_build (str): Target genome build for output
    - chroms (list[str]): List of chromosomes to write output for
    - coordinate: coordinate system string
    """

    vcfs = glob.glob(vcf)
    normalized_chroms = [normalize_chrom(c) for c in chroms]

    # Read sumstats with pandas; normalise chromosome column
    chr_remap = {"23": "X", "24": "Y", "25": "MT"}
    ssdf = pd.read_csv(
        ss, sep="\t", dtype=str,
        na_values=["NA", "NaN", "", "nan", "#NA"],
        low_memory=False,
    )
    if CHR_DSET in ssdf.columns:
        ssdf[CHR_DSET] = ssdf[CHR_DSET].astype(str).str.upper().map(
            lambda x: chr_remap.get(x, x)
        )
        ssdf = ssdf[ssdf[CHR_DSET].isin(normalized_chroms)]

    # handle the empty input file — still create output files so Nextflow output checks pass
    if ssdf.empty:
        print("No records in input summary statistics file.")
        for chrom in normalized_chroms:
            chrom_str = str(chrom).split(".")[0]
            open(f"{chrom_str}.merged", "w").close()
        open("unmapped", "w").close()
        return

    add_fields_if_missing(df=ssdf)
    rsid_mask = ssdf[RSID].str.startswith("rs").fillna(False)
    ssdf_with_rsid = ssdf[rsid_mask].copy()
    ssdf_without_rsid = ssdf[~rsid_mask].copy()
    header = list(ssdf.columns.values)

    print("starting rsid mapping")
    print("ssdf with rsid empty?: {}".format(ssdf_with_rsid.empty))

    # Use PyArrow dataset with is_in predicate pushdown — reads only rows whose
    # ID matches one of the rsIDs we need, never loading the full parquet into RAM.
    merged_vcf = pd.DataFrame()
    if not ssdf_with_rsid.empty:
        rsid_list = ssdf_with_rsid[RSID].dropna().tolist()

        ref_parts = []
        for vcf_path in vcfs:
            ref_table = ds.dataset(vcf_path).to_table(
                filter=ds.field("ID").isin(rsid_list),
                columns=["ID", "CHR", "POS"],
            )
            if ref_table.num_rows > 0:
                ref_parts.append(ref_table.to_pandas())

        if ref_parts:
            ref_df = pd.concat(ref_parts, ignore_index=True)
            # Deduplicate multi-chrom rsid hits — keep first (same as QUALIFY ROW_NUMBER()=1)
            ref_df = ref_df.drop_duplicates(subset=["ID"], keep="first")

            matched = ssdf_with_rsid.merge(ref_df, left_on=RSID, right_on="ID", how="inner")

            if not matched.empty:
                matched[CHR_DSET] = matched["CHR"].astype(str).str.replace(r"\..*$", "", regex=True)
                matched[BP_DSET] = matched["POS"].astype(str).str.replace(r"\..*$", "", regex=True)
                matched[HM_CC_DSET] = "rs"
                merged_vcf = matched[header + [HM_CC_DSET]]
                mapped_rsids = set(matched[RSID])
                ssdf_with_rsid = ssdf_with_rsid[~ssdf_with_rsid[RSID].isin(mapped_rsids)].copy()

    print("finished rsid mapping")
    # liftover the snps without rsids and those with unrecognised rsids 
    print("liftover remaining variants")
    ssdf = pd.concat([ssdf_with_rsid, ssdf_without_rsid])
    ssdf[HM_CC_DSET] = "lo"

    build_map = None
    if from_build != to_build:
        build_map = lft.LiftOver(lft.ucsc_release.get(from_build), lft.ucsc_release.get(to_build))
        
    @lru_cache(maxsize=100000)
    def cached_liftover(chrom, bp):
        try:
            return lft.map_bp_to_build_via_liftover(chromosome=chrom, bp=bp, build_map=build_map, coordinate=coordinate[0])
        except:
            return None
    # liftover the unmapped variants by chr and position
    if build_map:
        ssdf[BP_DSET] = [
            cached_liftover(str(chrom), str(bp)) if pd.notnull(chrom) and pd.notnull(bp) else None
            for chrom, bp in zip(ssdf[CHR_DSET], ssdf[BP_DSET])
            ]
    print("liftover complete")
    # merge "rs" and "lo" result to write the output
    combined_df = pd.concat([merged_vcf, ssdf], ignore_index=True)
    combined_df[CHR_DSET] = combined_df[CHR_DSET].astype("str").str.replace("\..*$","",regex=True)
    combined_df[BP_DSET] = combined_df[BP_DSET].astype("str").str.replace("\..*$","",regex=True)
    
    # 1. Write variants missing CHR or BP to "unmapped"
    unmapped_df = combined_df[combined_df[CHR_DSET].isnull() | combined_df[BP_DSET].isnull()].copy()
    unmapped_outfile = os.path.join("unmapped")
    unmapped_df.to_csv(unmapped_outfile, sep="\t", index=False, na_rep="NA")
    
    # 2. Write valid variants per chromosome, sorted by position so downstream
    #    tabix fetches are sequential (critical for performance on spinning disks).
    valid_df = combined_df.dropna(subset=[CHR_DSET, BP_DSET]).copy()
    valid_df["_bp_int"] = pd.to_numeric(valid_df[BP_DSET], errors="coerce")
    valid_df = valid_df.sort_values("_bp_int").drop(columns=["_bp_int"])

    for chrom in normalized_chroms:
        chrom_str = str(chrom).split(".")[0]
        out_path = os.path.join("{}.merged".format(chrom_str))
        chrom_df = valid_df[valid_df[CHR_DSET] == chrom]
        if not chrom_df.empty:
            chrom_df.to_csv(out_path, sep="\t", index=False, na_rep="NA")
        else:
            valid_df.head(0).to_csv(out_path, sep="\t", index=False)

def listify_string(string):
    """
    listify the input. If it's a list leave it.
    If it looks like a list literal (contains '[' and ']'), parse it.
    Otherwise treat the whole string as a single element.
    :param string:
    :return: a list
    """
    if type(string) is list:
        return string
    if type(string) is str:
        if "[" in string and "]" in string:
            new = string.replace(" ", "").replace("[", '["').replace("]", '"]').replace(",", '","')
            return literal_eval(new)
        return [string]
    return [str(string)]

def add_fields_if_missing(df):
    add_column_to_df(df=df, column=RSID)
    add_column_to_df(df=df, column=CHR_DSET)
    add_column_to_df(df=df, column=BP_DSET)

def add_column_to_df(df, column, value='NA'):
    if column not in df.columns:
        df[column] = value



def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-f', help='The name of the file to be processed', required=True)
    argparser.add_argument('-vcf', help='The name of the vcf file', required=True)
    argparser.add_argument('--log', help='The name of the log file')
    argparser.add_argument('-from_build', help='The original build e.g. "36" for NCBI36 or hg18', required=True)
    argparser.add_argument('-to_build', help='The latest (desired) build e.g. "38"', required=True)
    argparser.add_argument('-chroms', help='A chromosome or list of chromosomes to process', default=DEFAULT_CHROMS)
    argparser.add_argument('-coordinate', help='index', nargs='?', const="1-based", required=True)
    args = argparser.parse_args()

    ss = args.f
    vcf = args.vcf
    from_build = args.from_build
    to_build = args.to_build
    chroms = listify_string(args.chroms)
    coordinate = args.coordinate

    merge_ss_vcf(ss, vcf, from_build, to_build, chroms, coordinate)



if __name__ == "__main__":
    main()