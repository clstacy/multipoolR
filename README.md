# multipoolR: Multipool QTL Mapping in R

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**multipoolR** is an R package designed for Quantitative Trait Locus (QTL) mapping using data from bulk segregant analysis (BSA) experiments involving two pools sequenced with next-generation sequencing (NGS). It implements a statistical approach based on Kalman filtering and smoothing to estimate allele frequencies along the genome and calculates Logarithm of Odds (LOD) scores to identify potential QTLs.

The package is adapted from the concepts and methods described in Edwards & Gifford (2012) and the original Python implementation.

## Core Features

-   **Kalman Filter & Smoother:** Estimates underlying allele frequency trajectories and posterior estimates, accounting for recombination and sampling noise.
-   **Two Analysis Modes:**
    -   `replicates`: For comparing replicate pools selected under the *same* condition against the null hypothesis of 0.5 allele frequency.
    -   `contrast`: For comparing pools selected under *different* conditions, testing for significant allele frequency differences between them.
-   **LOD Score Calculation:** Identifies genomic regions associated with the trait difference based on likelihood ratios.
-   **Permutation Testing & FDR:** Empirically estimates significance (p-values) and controls the False Discovery Rate (q-values) using genome-wide permutations.
-   **Gene Annotation:** Automatically annotates significant QTL regions using standard Bioconductor annotation databases (e.g., for *Saccharomyces cerevisiae*).
-   **Genome-wide Visualization:** Generates publication-quality plots of allele frequencies and LOD scores using `ggplot2`.

## Installation

You can install the development version of `multipoolR` from GitHub using the `devtools` package:

``` r
# install.packages("devtools") # If you don't have devtools installed
devtools::install_github("clstacy/multipoolR")
```

## Dependencies

`multipoolR` requires R (\>= 4.0 recommended) and several packages from CRAN and Bioconductor.

-   **CRAN:**
    -   `ggplot2`
    -   `ggrepel`
    -   `checkmate`
    -   `lifecycle`
    -   `matrixStats`
    -   `future.apply` (Optional, only needed if using `parallel = TRUE` for permutations)
-   **Bioconductor:**
    -   `GenomicRanges`
    -   `GenomicFeatures`
    -   `AnnotationDbi`
    -   `GenomeInfoDb`
    -   `BiocGenerics`
    -   `S4Vectors`
    -   Specific annotation packages (e.g., `TxDb.Scerevisiae.UCSC.sacCer3.sgdGene`, `org.Sc.sgd.db` for yeast) are needed for the annotation feature and will be suggested during runtime if not installed.

You can install Bioconductor packages using `BiocManager`:

``` r
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("GenomicRanges", "GenomicFeatures", "AnnotationDbi", "GenomeInfoDb", "BiocGenerics", "S4Vectors"))
# BiocManager::install("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene") # Example yeast annotation
# BiocManager::install("org.Sc.sgd.db") # Example yeast annotation
```

## Input Data Format

Input data for each pool can be provided as either:

1.  A file path to a whitespace-delimited text file.
2.  A data frame already loaded into R.

The data must contain the following columns (default names shown, can be customized using `col_*` arguments):

-   `chr`: Chromosome identifier (e.g., "chrI", "chrII", "1", "2"). Must be consistent between pools and match the annotation database style if used.
-   `pos`: Genomic position (numeric base pair coordinate).
-   `a`: Read count supporting allele 'A' (numeric integer).
-   `b`: Read count supporting allele 'B' (numeric integer).

A header line is expected by default (`header = TRUE`). If no header is present, set `header = FALSE` and ensure columns are in the order `chr` (optional), `pos`, `a`, `b`. If the `chr` column is missing (3-column format), `assume_chr` will be used.

## Basic Usage

The main function is `multipool()`. Here's a basic example using the included example dataset.

``` r
library(multipoolR)
library(ggplot2) # For displaying the plot

# --- 1. Load Example Data ---
# This dataset is included with the package
data(multipoolR_example_data) # Loads multipoolR_example_pool1 and multipoolR_example_pool2

# --- 2. Run multipool analysis ---
# Ensure required annotation packages are available if you want annotation
# library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene) # for yeast
# library(org.Sc.sgd.db) # for yeast

results_mp <- multipool(
  pool1 = multipoolR_example_pool1,    # Use the loaded data
  pool2 = multipoolR_example_pool2,    # Use the loaded data
  N = 50,              # Effective number of individuals (example value)
  mode = "contrast", # Test for difference in alleles between pools
  res = 1000,          # Bin size in base pairs (adjust as needed)
  cM = 3300,           # Estimated bp per centiMorgan (yeast default)
  nperm = 10,        # Number of permutations for FDR (use >= 1000 for publication)
  q_threshold = 0.05,  # Q-value threshold for significance
  plot = TRUE,         # Generate the genome-wide plot
  seed = 42            # Set seed for reproducible permutations
  # txdb = TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, # Optional: Provide annotation DBs
  # orgdb = org.Sc.sgd.db
)

# --- 3. Explore Results ---

# View the results table (includes LOD, p-values, q-values)
print(head(results_mp$results))

# Check the LOD score threshold corresponding to the q-value cutoff
if (!is.null(results_mp$fdr_lod_threshold) && is.finite(results_mp$fdr_lod_threshold)) {
    print(paste("LOD threshold for q <", results_mp$parameters$q_threshold, ":",
                round(results_mp$fdr_lod_threshold, 3)))
} else {
    print("No significant LOD threshold found at the specified q-value.")
}


# Display the plot (if plot = TRUE)
# The plot shows observed (points) and fitted (lines) frequencies for pool1 (blue)
# and pool2 (red), along with the LOD score (black line, secondary axis).
# A dashed orange line indicates the FDR threshold. Significant genes may be labeled.
if (!is.null(results_mp$plot)) {
  print(results_mp$plot)
}

# View annotated genes in significant regions (if annotation DBs provided)
if (!is.null(results_mp$annotated_genes) && nrow(results_mp$annotated_genes) > 0) {
  print(head(results_mp$annotated_genes))
} else {
    print("No annotated genes found or annotation was skipped.")
}

# Access parameters used
print(results_mp$parameters)
```

## Key Parameters

-   `pool1`, `pool2`: Input data (file paths or data frames).
-   `N`: Effective population size (crucial for variance estimation).
-   `mode`: `"replicates"` or `"contrast"`.
-   `res`: Bin resolution (bp). Smaller bins give higher resolution but may have lower counts.
-   `cM`: Genetic map density (bp/cM). Affects the recombination rate `p` used in the Kalman filter.
-   `nperm`: Number of permutations. **Set \> 0** (e.g., 1000+) for reliable p-values **and q-values.** `nperm = 0` runs the analysis faster but provides no statistical significance assessment beyond raw LOD scores.
-   `q_threshold`: Significance cutoff for FDR (q-value). Used for determining the `fdr_lod_threshold` and selecting genes for annotation/plotting.
-   `filter`: Apply default filters for low/high coverage and fixated markers (recommended).
-   `txdb`, `orgdb`: Bioconductor annotation database objects (required for annotation).
-   `col_chr`, `col_pos`, `col_a`, `col_b`: Customize input column names.
-   `header`, `assume_chr`: Control file parsing behavior.
-   `parallel`, `seed`: Control permutation execution.
-   `verbose`: Logical, control whether progress messages are printed (default: TRUE).

## Citation

If you use `multipoolR` in your research, please cite the paper describing the statistical method:

Edwards, M. D., & Gifford, D. K. (2012). High-resolution genetic mapping with pooled sequencing. *BMC bioinformatics*, 13(Suppl 6), S8. <https://doi.org/10.1186/1471-2105-13-S6-S8>

## License

This package is licensed under the MIT License. See the LICENSE file for details.
