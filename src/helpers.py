#!/usr/bin/env python3

"""
Helper function for analyzing Brian Peterson's data
"""

import glob
import pandas as pd
import numpy as np
from statsmodels.stats import rates

EPSILON = np.finfo(float).eps


def get_array(info_dict, aa_num):

    length = 0
    for region,reg_info in info_dict.items():
        length += reg_info["ratio"].shape[1]

    ratios = []
    rat_var = []
    pvals = []
    p_var = []

    for region,reg_info in info_dict.items():
        ratios.append(reg_info["ratio"])
        rat_var.append(reg_info["rat_var"])
        pvals.append(reg_info["pval"])
        p_var.append(reg_info["pval_var"])

    rat_arr = np.concatenate(ratios, axis=1)
    rat_var_arr = np.concatenate(rat_var, axis=1)
    pval_arr = np.concatenate(pvals, axis=1)
    pval_var_arr = np.concatenate(p_var, axis=1)

    return (rat_arr, rat_var_arr, pval_arr, pval_var_arr)


def calculate_wt(group):
    exposure = group["exposure"][0]
    var_count = group["count"].sum()
    return exposure - var_count


def filter_files(file_list, term):
    return [_ for _ in file_list if term in _]


def read_data_to_dict(file_list, regions, reps, lut, name_col, aa_row_lut):

    df_dict = {}

    for region in regions:

        start,stop = lut[region]
        npos = stop - start + 1

        if not region in df_dict:
            df_dict[region] = {}
 
        wt_aa_list = []
        rates_list = []
        counts_list = []
        exposure_list = []
        for rep in reps:
            rep_files = filter_files(file_list, rep)
            X = np.zeros((len(aa_row_lut), npos))
            this_file = filter_files(rep_files, region)
            if len(this_file) == 0:
                print(f"WARNING: no file found for {region}, {rep}. Skipping it.")
                continue
            if len(this_file) > 1:
                raise(Exception(
                    f"Expecting one file, but multiple in list after filtering "\
                    f"for region {region} and replicate {rep}.\n"\
                    f"File list:\n{this_file}"
                ))
            print(f"Parsing mutation counts in {this_file}")
            df = pd.read_table(this_file[0], header=0)
            df["replicate"] = rep
            df["region"] = region
            df.rename(
                columns = {
                    "Amino Acid Position-Start=1":"position",
                    "1-mutation_read_count":"singles",
                    "WT_read_count":"wt_count",
                    "Ref_AA":"ref_aa",
                },
                inplace = True,
            )
            df["exposure"] = df["singles"] + df["wt_count"]
            df = df.loc[
                (df["position"] >= start)
                & (df["position"] <= stop)
            ,:]
            df["outcome"] = name_col

            df_long = pd.melt(
                df,
                id_vars = [
                    "position",
                    "outcome",
                    "ref_aa",
                    #"exposure",
                    "replicate",
                    "region",
                ],
                value_vars = [
                    "A", "G", "I", "L", "P", "V", "F", "W", "Y",
                    "D", "E", "R", "H",
                    "K", "S", "T", "C", "M", "N", "Q", "*", "X",
                    # "-",
                ],
                value_name = "count",
                var_name = "var_aa",
            )

            for i,row in df_long.iterrows():
                x_row_idx = aa_row_lut[row.var_aa]
                x_col_idx = row.position - start
                X[x_row_idx, x_col_idx] = row["count"]

            rate = X / df["exposure"].values[None]
            rates_list.append(rate)
            counts_list.append(X)
            exposure_list.append(df["exposure"].values[None])

            #mut_df = df_long[
            #    df_long["ref_aa"] != df_long["var_aa"]
            #].sort_values(["position","var_aa"])

            #position_list = []
            #var_aa_list = []
            #interaction_list = []
            #for i,row in df_long.iterrows():
            #    for j in range(row["count"]):
            #        position = row.position
            #        if row.ref_aa == row.var_aa:
            #            position = 0
            #        position_list.append(position)
            #        var_aa_list.append(row.var_aa)
            #        interaction_list.append(f"{position}_{row.var_aa}")

            #read_df = pd.DataFrame({
            #    "position":position_list,
            #    "var":var_aa_list,
            #    "interaction":interaction_list,
            #})
            #X = pd.get_dummies(
            #    read_df,
            #    columns=["position","var","interaction"],
            #    sparse=True,
            #)
            #X["intercept"] = 1

            #df_dict[region][rep] = {
            #    "X":X,
            #    "exposure":df["exposure"].values[None,:],
            #    "rate": X / df["exposure"].values[None],
            #}

        df_dict[region] = {
            "rates": np.stack(rates_list, axis=-1),
            "counts": np.stack(counts_list, axis=-1),
            "exposures": np.stack(exposure_list, axis=-1),
            "ref_aa": df["ref_aa"].values,
        }
    

    return df_dict


def set_up_dataframe(df_dict):

    first = True
    for rep,rep_vals in df_dict.items():
        for reg,reg_df in rep_vals.items():
            if first:
                df = reg_df.copy()
                first = False
            else:
                df = pd.concat([df, reg_df])

    return df_long


def calc_enrichment(count_a, exposure_a, count_base, exposure_base, direction):
    ratios = np.zeros_like(count_a)
    stats = np.zeros_like(count_a)
    pvals = np.zeros_like(count_a)
    for i in range(count_a.shape[0]):
        for j in range(count_a.shape[1]):
            result = rates.test_poisson_2indep(
                count_a[i,j],
                exposure_a[0,j],
                count_base[i,j],
                exposure_base[0,j],
                alternative = direction,
            )
            ratios[i,j] = result.ratio
            stats[i,j] = result.statistic
            pvals[i,j] = result.pvalue
    return ( ratios, stats, pvals )


def jack_var(theta_hat, thetas, alphas):
    R = len(alphas)
    var = np.sum([(alphas[r] * (thetas[r] - theta_hat)**2) for r in range(R)])
    return var


def jackknife(sel_counts, sel_exp, input_counts, input_exp, direction):
    
    all_counts = np.concatenate((sel_counts, input_counts), axis=2)
    all_exps = np.concatenate((sel_exp, input_exp), axis=2)
    sel_reps = sel_counts.shape[2]
    in_reps = input_counts.shape[2]
    strata = np.array(
        [0 for _ in range(sel_reps)] + [1 for _ in range(in_reps)]
    )
    nh = {0: sel_reps, 1: in_reps}
    jack_rep_num = all_counts.shape[2]
    W = [1/sel_reps for _ in range(sel_reps)] + [1/in_reps for _ in range(in_reps)]
    jack_ratios = []
    jack_pvals = []
    alphas = []

    for jack_idx in range(jack_rep_num):

        jack_weights = np.zeros(jack_rep_num)
        jack_stratum = strata[jack_idx]
        alpha = (nh[jack_stratum] - 1) / nh[jack_stratum]
        alphas.append(alpha)

        # set jackknife weights for this jack replicate
        for w_i,w in enumerate(jack_weights):
            stratum = strata[w_i]
            if w_i == jack_idx:
                continue
            if stratum == jack_stratum:
                jack_weights[w_i] = W[w_i]/alpha
            else:
                jack_weights[w_i] = W[w_i]

        jack_counts = all_counts * jack_weights[None,None,:]
        sel_jack_count = jack_counts[:,:,strata==0].sum(axis=2)
        in_jack_count = jack_counts[:,:,strata==1].sum(axis=2)
        jack_exp = all_exps * jack_weights[None,None,:]
        sel_jack_exp = jack_exp[:,:,strata==0].sum(axis=2)
        in_jack_exp = jack_exp[:,:,strata==1].sum(axis=2)
        ratio,stat,pval = calc_enrichment(
            sel_jack_count,
            sel_jack_exp,
            in_jack_count,
            in_jack_exp,
            direction,
        )
        jack_ratios.append(ratio)
        jack_pvals.append(pval)

    jack_rat_arr = np.stack(jack_ratios, axis=-1)
    jack_pval_arr = np.stack(jack_pvals, axis=-1)

    enrichment = np.mean(jack_rat_arr, axis=2)
    pval = np.mean(jack_pval_arr, axis=2)
    enrich_var = np.zeros_like(enrichment)
    pval_var = np.zeros_like(pval)
    for i in range(enrich_var.shape[0]):
        for j in range(enrich_var.shape[1]):
            enrich_var[i,j] = jack_var(
                enrichment[i,j],
                jack_rat_arr[i,j,:],
                alphas,
            )
            pval_var[i,j] = jack_var(
                pval[i,j],
                jack_pval_arr[i,j,:],
                alphas,
            )

    result_dict = {
        "ratio": enrichment,
        "rat_var": enrich_var,
        "pval": pval,
        "pval_var": pval_var,
    }

    return result_dict


def get_enrichments(sel_data, input_data, direction="two-sided"):

    enrich_dict = {}

    for reg,sel_info in sel_data.items():

        ref_aa = sel_info["ref_aa"]
        sel_counts = sel_info["counts"]
        input_counts = input_data[reg]["counts"]
        sel_exposures = sel_info["exposures"]
        input_exposures = input_data[reg]["exposures"]

        enrichments = jackknife(
            sel_counts,
            sel_exposures,
            input_counts,
            input_exposures,
            direction=direction,
        )
        enrichments["ref_aa"] = ref_aa
        enrich_dict[reg] = enrichments

    return enrich_dict

