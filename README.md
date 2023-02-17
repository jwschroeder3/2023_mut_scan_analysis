# Hrd1 variant analysis

Code to perform analysis of Hrd1 mutational scanning data is in `src`.
To set up the conda environment required to run the code in `src/driver.py`,
enter `src`, then, assuming you have `conda` working on your system,
run `conda env create -f conda_environment.yaml`.

This will set up the conda environment `mutscan`.

## Analysis

Enter the top directory of this repository in a bash terminal
and run the following to reproduce our analysis:

```bash
conda activate mutscan
python src/driver.py --in-direc counts --out-direc output
```

For each Hrd1 phenotype selected (ldead, mdead, wt), the above code performs
many one-tailed rate-ratio tests, one test for every combination
of amino acid residue and position in the given category
of Hrd1 phenotype. Jackknife sampling was performed to incorporate
replicate information into our estimates of p-values, log(fold-enrichment),
and variance in our point estimates.

## Expected results

In the `output` directory you should see several csv files with
the results of having run the above code.

The prefix `ldead`, `mdead` or `wt` indicates the a given file
includes mutation enrichment results for the indicated category
of Hrd1 phenotype described in the accompanying publication.

Below we describe the contents of files with each indicated suffix:

* `*_BH.csv` - This is a primary file of interest, as it contains adjusted p-values to account for multiple hypothesis testing. We used the Benjamini-Hochberg method to adjust p-values.
* `*_pvals.csv` - jackknife mean p-values for one-tailed rate-ratio test testing whether mutation to a gien amino acid at a given Hrd1 position is enriched in the given Hrd1 phenotype (ldead, mdead, wt)
* `*_pvals_variance.csv` - jackknife variance in the p-value estimate
* `*_ratios.csv` - jackknife mean ratio for enrichment of mutations to a given amino acid at a given position of Hrd1 in the given Hrd1 phenotype (note that any ratio less than 1.0 was clipped to 1.0, as we were only interested in enrichments for this analysis, as opposed to depletions)
* `*_ratios_log.csv` - log2-transformed ratios
* `*_ratio_variance` - jackknife variance for the enrichment estimate

