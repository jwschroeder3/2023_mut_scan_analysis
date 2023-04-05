#!/usr/bin/env python3

"""
Driver script for doing rate ratio test on Brian Peterson's
Hrd1 data.
"""

import os
import sys
import numpy as np
import pandas as pd
import argparse
import glob
from statsmodels.stats.multitest import multipletests

import helpers as h

def write_to_file(fname, arr, header, ref_aa):
    with open(fname, "w") as outf:
        outf.write(header)
        for c_idx in range(arr.shape[1]):
            data = ','.join([str(_) for _ in arr[:,c_idx]])
            outf.write(f"{data},{c_idx+1},{ref_aa[c_idx]}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-direc', action='store',
        type=str, required=True,
        help='path to directory containing input data')
    parser.add_argument('--out-direc', action='store',
        type=str, required=True,
        help='path to directory containing outputs')
    #parser.add_argument('--read-class', action='store',
    #    type=str, required=True,
    #    help=f"Controls whether to use all-, single-, or multi-"\
    #        f"mutation reads for number of mutations")

    args = parser.parse_args()

    in_direc = args.in_direc
    out_direc = args.out_direc
    #read_class = args.read_class
    read_class = "Single"

    reps = [f"Rep{_+1}" for _ in range(4)]
    regions = [f"Region{_+1}" for _ in range(5)]
    region_key = {
        "Region1": (1,   110),
        "Region2": (111, 220),
        "Region3": (221, 330),
        "Region4": (331, 440),
        "Region5": (441, 551),
    }
    phenotypes = ["LDead", "MDead", "WT"]
    amino_acids = [
        "A","G","I","L","P","V","F","W","Y","D","E",
        "R","H","K","S","T","C","M","N","Q","*","X",
    ]

    row_aa_lut = {aa:i for i,aa in enumerate(amino_acids)}

    all_files = glob.glob(os.path.join(in_direc, f"*_{read_class}.tsv"))


    ldead_files = h.filter_files(all_files, phenotypes[0])
    wt_files = h.filter_files(all_files, phenotypes[2])
    input_files = h.filter_files(all_files, "Input")


    ldead_data = h.read_data_to_dict(
        ldead_files,
        regions,
        reps,
        region_key,
        "ldead",
        row_aa_lut,
    )
    wt_data = h.read_data_to_dict(
        wt_files,
        regions,
        reps,
        region_key,
        "wt",
        row_aa_lut,
    )
    input_data = h.read_data_to_dict(
        input_files,
        regions,
        reps,
        region_key,
        "input",
        row_aa_lut,
    )

    ref_aa = []
    for reg,reg_info in input_data.items():
        ref_aa.extend(list(reg_info["ref_aa"]))

    ldead_enrichment = h.get_enrichments(
        ldead_data,
        input_data,
        "larger",
    )
    wt_enrichment = h.get_enrichments(
        wt_data,
        input_data,
        "two-sided",
    )

    aa_row_lut = {row_idx:aa for aa,row_idx in row_aa_lut.items()}

    ldead_rats,ldead_rat_var,ldead_pvals,ldead_p_var = h.get_array(
        ldead_enrichment,
        len(aa_row_lut),
    )
    wt_rats,wt_rat_var,wt_pvals,wt_p_var = h.get_array(
        wt_enrichment,
        len(aa_row_lut),
    )

    arr_shape = wt_rats.shape

    header = f"{','.join(amino_acids)},position,ref_aa\n"

    out_fname = f"output/wt_Single_ratios.csv"
    write_to_file(out_fname, wt_rats, header, ref_aa)
    out_fname = f"output/wt_Single_ratios_log.csv"
    wt_log_rats = np.log2(wt_rats)
    write_to_file(out_fname, wt_log_rats, header, ref_aa)
    out_fname = f"output/wt_Single_ratio_variance.csv"
    write_to_file(out_fname, wt_rat_var, header, ref_aa)
    out_fname = f"output/wt_Single_pvals.csv"
    write_to_file(out_fname, wt_pvals, header, ref_aa)
    out_fname = f"output/wt_Single_pval_variance.csv"
    write_to_file(out_fname, wt_p_var, header, ref_aa)
    p_vals = wt_pvals.flatten()
    p_vals[np.isnan(p_vals)] = 1
    wt_bh = multipletests(p_vals, method="fdr_bh")[1]
    wt_q = np.reshape(wt_bh, arr_shape)
    out_fname = f"output/wt_Single_BH.csv"
    write_to_file(out_fname, wt_q, header, ref_aa)
 
 
    out_fname = f"output/ldead_Single_ratios.csv"
    write_to_file(out_fname, ldead_rats, header, ref_aa)
    out_fname = f"output/ldead_Single_ratios_log.csv"
    ldead_rats[ldead_rats < 1.0] = 1.0
    ldead_log_rats = np.log2(ldead_rats)
    write_to_file(out_fname, ldead_log_rats, header, ref_aa)
    out_fname = f"output/ldead_Single_ratio_variance.csv"
    write_to_file(out_fname, ldead_rat_var, header, ref_aa)
    out_fname = f"output/ldead_Single_pvals.csv"
    write_to_file(out_fname, ldead_pvals, header, ref_aa)
    out_fname = f"output/ldead_Single_pval_variance.csv"
    write_to_file(out_fname, ldead_p_var, header, ref_aa)
    p_vals = ldead_pvals.flatten()
    p_vals[np.isnan(p_vals)] = 1
    ldead_bh = multipletests(p_vals, method="fdr_bh")[1]
    ldead_q = np.reshape(ldead_bh, arr_shape)
    out_fname = f"output/ldead_Single_BH.csv"
    write_to_file(out_fname, ldead_q, header, ref_aa)
 


if __name__ == '__main__':
    main()


