#' multipoolR: Multipool QTL Mapping in R
#'
#' A comprehensive package to perform multipool QTL scans using a Kalman filter
#' approach, permutation testing for significance (calculating q-values via FDR),
#' yeast gene annotation for significant QTLs, and genome-wide visualization of results.
#' Adapted from the original Python implementation by Edward and Giffords.
#'
#' @docType _PACKAGE
#' @name multipoolR
#'
#' @import ggplot2
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import ggrepel
#' @import checkmate
#' @importFrom AnnotationDbi select keytypes
#' @importFrom lifecycle deprecate_warn
#' @importFrom matrixStats colMaxs
#' @importFrom future.apply future_lapply
#' @importFrom stats dnorm median quantile qnorm setNames aggregate na.omit optimize p.adjust
#' @importFrom utils read.table write.table head tail
#' @importFrom graphics abline hist legend par plot text
#' @importFrom grDevices dev.off pdf png
#' @importFrom GenomeInfoDb mapSeqlevels seqlevelsStyle seqlevels seqlevels<- seqlevelsStyle<-
#' @importFrom S4Vectors mcols subjectHits queryHits
#' @importFrom BiocGenerics start end order
#'
#' @section Main Function:
#' \code{\link{multipool}}
#'
#' @section Data Format:
#' Input data (either files or data frames) should contain columns:
#' \describe{
#'   \item{chr}{Chromosome identifier (e.g., "chrI", "chrII"). Must be consistent between pools.}
#'   \item{pos}{Genomic position (numeric).}
#'   \item{a}{Read count supporting allele 'A' (numeric).}
#'   \item{b}{Read count supporting allele 'B' (numeric).}
#' }
#'
#' @section Algorithm Details:
#' The core algorithm uses a Kalman filter and smoother to estimate the underlying
#' allele frequency trajectory along the genome for each pool, accounting for
#' recombination and sampling variance. LOD scores are calculated based on the
#' likelihood ratio comparing a model with a QTL at a specific location versus
#' a null model (no QTL), maximized over possible true allele frequencies using numerical optimization.
#' Genome-wide permutation tests can be performed to estimate empirical p-values and q-values (False Discovery Rate).
#'
#' @references
#' Magwene, P. M., Willis, J. H., & Kelly, J. K. (2011). The statistics of bulk segregant analysis using next generation sequencing. PLoS computational biology, 7(11), e1002255.
NULL
# Needed to avoid R CMD check NOTES about undefined global variables
utils::globalVariables(c(
  "pos", "freq1_obs", "freq2_obs", "LOD_scaled", "chr", "freq", "pool",
  "ylo", "yhi", "alpha_tag", "GENENAME", "start", "lod_scaled", "x",
  "peak_bin_pos", "yintercept_scaled", "threshold_name", "pos_mb", # Added vars used in ggplot aes
  "bin_start", # Added var used in aggregate
  "gene_id", # Added var used in annotation split/merge
  "peak_bin_lod", # Added var used in plotting filter
  "threshold_value", # Added for plotting thresholds
  "q_value",
  "min_q_value" # Added for annotation/plotting
))

# ---- Constants ----
DELTA_X <- 0.0025 # Step size for integrating allele frequencies in LOD calc (used for range)
N_GENES_LABEL <- 3 # Max number of genes to label per chromosome peak region
# Default q-value threshold for significance when annotating/plotting if nperm > 0
DEFAULT_Q_THRESHOLD <- 0.05

msg <- function(..., level = 1) {
  ## look for a package-wide option; falls back to TRUE
  if (isTRUE(getOption("multipool.verbose", TRUE)) && level >= 1)
    base::message(...)
}

# ---- Main Multipool Function ----

#' Perform Multipool Genome Scan
#'
#' This is the main function to run the multipool analysis. It loads data for
#' two pools, performs QTL mapping using a Kalman filter approach, optionally
#' runs permutation tests to estimate significance (p-values and q-values),
#' annotates significant QTL regions with yeast gene information, and generates
#' a genome-wide plot of the results.
#'
#' @param pool1 A file path (character string) to the first pool's data, or a
#'   data frame containing the required columns (`chr`, `pos`, `a`, `b`).
#' @param pool2 A file path (character string) to the second pool's data, or a
#'   data frame containing the required columns (`chr`, `pos`, `a`, `b`).
#' @param N Numeric. The effective number of individuals contributing to each
#'   sequencing pool. This accounts for the variance introduced during pooling
#'   and sequencing.
#' @param mode Character string. Specifies the analysis mode:
#'   \describe{
#'     \item{"replicates"}{Assumes the two pools are biological or technical
#'       replicates selected under the *same* condition. The LOD score tests
#'       for deviations from the expected 0.5 allele frequency in *both* pools
#'       simultaneously.}
#'     \item{"contrast"}{Assumes the two pools were selected under *different*
#'       conditions. The LOD score tests for significant *differences* in
#'       allele frequency between the two pools.}
#'   }
#' @param res Numeric. The desired bin resolution in base pairs for analysis
#'   (default: 100). Marker data within each bin is aggregated.
#' @param cM Numeric. An estimate of the average number of base pairs per
#'   centiMorgan for the organism (default: 3300 for yeast). Used to calculate
#'   the recombination rate between adjacent bins.
#' @param filter Logical. If TRUE (default), apply filtering to remove markers
#'   with zero counts for either allele ('fixated' markers) and bins with
#'   extremely high coverage (potential PCR duplicates/mapping artifacts).
#' @param nperm Integer. The number of permutations to perform for establishing
#'   significance (p-values and q-values) (default: 0, no permutations). Setting
#'   this to > 0 (e.g., 1000) is recommended for robust analysis.
#' @param alpha Deprecated. Significance level is now typically determined by choosing a
#'   q-value threshold (e.g., q < 0.05) on the results.
#' @param plot Logical. If TRUE (default), generate a genome-wide plot of the
#'   results using ggplot2.
#' @param txdb A TxDb object for gene annotations (default:
#'   TxDb.Scerevisiae.UCSC.sacCer3.sgdGene). See GenomicFeatures package.
#' @param orgdb An OrgDb object for gene annotations (default: org.Sc.sgd.db).
#'   See AnnotationDbi package.
#' @param q_threshold Numeric. The q-value threshold for calling significance
#'   in annotation and highlighting genes on the plot when `nperm > 0`.
#'   (Default: 0.05).
#'
#' @return A list containing:
#'   \describe{
#'     \item{results}{A data frame with genome-wide results, including columns
#'       for chromosome (`chr`), bin start position (`pos`), observed allele
#'       frequencies (`freq1_obs`, `freq2_obs`), smoothed posterior allele
#'       frequencies (`freq1_fit`, `freq2_fit`), LOD scores (`LOD`), MLE
#'       allele frequencies or differences (`mu_MLE_adj`), and if `nperm > 0`,
#'       empirical p-values (`p_value`) and q-values (`q_value`).}
#'     \item{fdr_lod_threshold}{The LOD score threshold corresponding to the chosen `q_threshold`,
#'       if permutations were run and significant results found. Otherwise NULL.}
#'     \item{plot}{A ggplot object containing the genome-wide plot. NULL if
#'       `plot` was FALSE.}
#'     \item{annotated_genes}{A data frame containing information about genes
#'       overlapping bins. If `nperm > 0`, this includes genes overlapping bins
#'       meeting the `q_threshold`. If `nperm = 0`, this includes
#'       genes overlapping *any* bin. NULL if annotation databases unavailable.}
#'     \item{parameters}{A list storing the key parameters used for the analysis.}
#'   }
#' @export
#' @examples
#' \dontrun{
#' # --- Example with dummy data ---
#'
#' # Create dummy data frames
#' make_dummy_pool <- function(n_markers = 5000, chr = "chrI") {
#'   pos <- sort(sample(1:200000, n_markers, replace = TRUE))
#'   a <- rpois(n_markers, lambda = 20)
#'   b <- rpois(n_markers, lambda = 20)
#'   # Introduce a fake QTL
#'   qtl_region <- pos > 80000 & pos < 120000
#'   a[qtl_region] <- rpois(sum(qtl_region), lambda = 35)
#'   b[qtl_region] <- rpois(sum(qtl_region), lambda = 10)
#'   data.frame(chr = chr, pos = pos, a = a, b = b)
#' }
#'
#' pool1_df <- make_dummy_pool(chr = "chrI")
#' pool2_df <- make_dummy_pool(chr = "chrI") # Replicate
#'
#' # Combine into two chromosomes for a more realistic example
#' pool1_combined <- rbind(pool1_df, make_dummy_pool(chr="chrII", n_markers=3000))
#' pool2_combined <- rbind(pool2_df, make_dummy_pool(chr="chrII", n_markers=3000))
#'
#' # Ensure required packages are loaded (needed for default annotation DBs)
#' # library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
#' # library(org.Sc.sgd.db)
#'
#' # Run analysis in "replicates" mode with permutations for FDR
#' results_rep <- multipool(
#'   pool1 = pool1_combined,
#'   pool2 = pool2_combined,
#'   N = 50,
#'   mode = "replicates",
#'   res = 500,
#'   cM = 3300,
#'   nperm = 1000, # Use >= 1000 for reliable FDR
#'   q_threshold = 0.10 # Use a 10% FDR threshold
#' )
#'
#' # View results with p/q-values
#' print(head(results_rep$results))
#' # Print the LOD threshold corresponding to the FDR cutoff
#' print(paste("LOD threshold for q <", results_rep$parameters$q_threshold, ":",
#'             round(results_rep$fdr_lod_threshold, 3)))
#' # View plot (will highlight genes with q < 0.1 and show threshold line)
#' if (!is.null(results_rep$plot)) {
#'   print(results_rep$plot)
#' }
#' print(results_rep$annotated_genes) # Genes meeting q-threshold
#'
#' # Run analysis without permutations (nperm = 0)
#' results_no_perm <- multipool(
#'   pool1 = pool1_combined,
#'   pool2 = pool2_combined,
#'   N = 50,
#'   mode = "replicates",
#'   res = 500,
#'   cM = 3300,
#'   nperm = 0
#' )
#'
#' # View plot (should label genes near peaks)
#' if (!is.null(results_no_perm$plot)) {
#'   print(results_no_perm$plot)
#' }
#' print(results_no_perm$annotated_genes) # All overlapping genes
#'
#' } # end dontrun
multipool <- function(pool1, pool2, N,
                      mode = c("replicates", "contrast"),
                      res = 100, cM = 3300, filter = TRUE,
                      nperm = 0, alpha = NULL, plot = TRUE, # alpha is deprecated
                      col_chr   = "chr",
                      col_pos   = "pos",
                      col_a     = "a",
                      col_b     = "b",
                      assume_chr = "chrI",
                      header    = TRUE,
                      seed = NULL, parallel = FALSE,
                      txdb = NULL, orgdb = NULL,
                      q_threshold = DEFAULT_Q_THRESHOLD,
                      verbose = TRUE) { # Added q_threshold



  # --- Input Validation and Setup ---
  mode <- match.arg(mode)

  # ---- rigorous argument checks ----
  checkmate::assert_number(N,    lower = 1, finite = TRUE)
  checkmate::assert_number(res,  lower = 1, finite = TRUE)
  checkmate::assert_number(cM,   lower = 1, finite = TRUE)
  checkmate::assert_int   (nperm, lower = 0)
  checkmate::assert_number(q_threshold, lower = 0, upper = 1)
  checkmate::assert_flag  (filter)
  checkmate::assert_flag  (plot)
  if (!is.null(alpha))
    lifecycle::deprecate_warn("0.2.0", "multipool(alpha)", "multipool(q_threshold)")


  # Load annotation databases if needed later
  loaded_txdb <- NULL
  loaded_orgdb <- NULL
  if (plot || nperm > 0) {
    if (is.null(txdb)) {
      if (requireNamespace("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", quietly = TRUE)) {
        loaded_txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene::TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
        msg("Using default yeast TxDb: TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
      } else {
        warning("Default TxDb package 'TxDb.Scerevisiae.UCSC.sacCer3.sgdGene' not found. Annotation will be skipped. Install it or provide a TxDb object to the 'txdb' argument.", call. = FALSE)
      }
    } else {
      loaded_txdb <- txdb
    }
    if (is.null(orgdb)) {
      if (requireNamespace("org.Sc.sgd.db", quietly = TRUE)) {
        loaded_orgdb <- org.Sc.sgd.db::org.Sc.sgd.db
        msg("Using default yeast OrgDb: org.Sc.sgd.db")
      } else {
        warning("Default OrgDb package 'org.Sc.sgd.db' not found. Annotation will be skipped. Install it or provide an OrgDb object to the 'orgdb' argument.", call. = FALSE)
      }
    } else {
      loaded_orgdb <- orgdb
    }
  }


  # Calculate recombination fraction between adjacent bins
  p_rec <- (res / 100) / cM
  if (p_rec >= 0.5) {
    warning(paste("Calculated recombination rate p =", round(p_rec, 3),
                  "is >= 0.5. Check 'res' and 'cM' values.",
                  "Setting p = 0.499 to proceed."), call. = FALSE)
    p_rec <- 0.499
  } else if (p_rec <= 0) {
    warning(paste("Calculated recombination rate p =", round(p_rec, 3),
                  "is <= 0. Check 'res' and 'cM' values.",
                  "Setting p = 1e-9 to proceed."), call. = FALSE)
    p_rec <- 1e-9
  }


  # --- Data Loading and Preprocessing ---
  msg("Loading and preprocessing data...")
  loaded_data <- internal_load_two_pools(pool1, pool2, res, filter,
                                         col_chr   = col_chr,
                                         col_pos   = col_pos,
                                         col_a     = col_a,
                                         col_b     = col_b,
                                         assume_chr = assume_chr,
                                         header    = header)

  if (length(loaded_data$chrom_data) == 0) {
    stop("No valid data loaded after processing. Check input files/data frames and parameters.")
  }

  # --- Core Computation (per chromosome) ---
  msg("Running core computation (Kalman filter and LOD calculation)...")
  all_results_list <- list()
  chrom_order <- loaded_data$common_chroms

  for (chrom in chrom_order) {
    if (!chrom %in% names(loaded_data$chrom_data)) next
    chr_data <- loaded_data$chrom_data[[chrom]]
    msg(paste("Processing chromosome:", chrom, "(", chr_data$T_, "bins )"))
    if (chr_data$T_ < 2) {
      warning(paste("Skipping chromosome", chrom, "due to insufficient bins (< 2) after processing."), call. = FALSE)
      next
    }
    comp_res <- internal_do_computation(
      y1 = chr_data$y1, y_var1 = chr_data$y_var1, d1 = chr_data$d1,
      y2 = chr_data$y2, y_var2 = chr_data$y_var2, d2 = chr_data$d2,
      T_ = chr_data$T_, N = N, p = p_rec, mode = mode
    )
    bin_starts <- chr_data$bins[-length(chr_data$bins)]
    results_df_chr <- data.frame(
      chr = chrom,
      pos = bin_starts,
      freq1_obs = chr_data$y1 / chr_data$d1,
      freq2_obs = chr_data$y2 / chr_data$d2,
      freq1_fit = comp_res$mu_pstr1 / N,
      freq2_fit = comp_res$mu_pstr2 / N,
      LOD = comp_res$LOD,
      mu_MLE_adj = comp_res$mu_MLE / N,
      stringsAsFactors = FALSE
    )
    results_df_chr$freq1_obs[chr_data$d1 == 0] <- NA
    results_df_chr$freq2_obs[chr_data$d2 == 0] <- NA
    all_results_list[[chrom]] <- results_df_chr
  }

  if (length(all_results_list) == 0) {
    stop("Computation failed for all chromosomes. Check input data quality and parameters.")
  }

  results_df_genome <- do.call(rbind, all_results_list[chrom_order])
  rownames(results_df_genome) <- NULL


  # --- Permutation Testing & FDR Calculation ---
  perm_LODs_all <- NULL # Store all permuted LODs
  fdr_lod_threshold <- NULL # Store the calculated LOD threshold for FDR
  if (nperm > 0) {
    msg(paste("Running genome-wide permutation test with nperm =", nperm, "..."))
    perm_LODs_all <- internal_multipool_permutation_fdr( # Call the FDR version
      chrom_data_list = loaded_data$chrom_data,
      N = N,
      p = p_rec,
      mode = mode,
      nperm = nperm,
      seed                 = seed,
      parallel             = parallel
    )

    if (!is.null(perm_LODs_all) && length(perm_LODs_all) > 0) {
      msg("Calculating p-values and q-values...")
      observed_LODs <- results_df_genome$LOD
      valid_obs_indices <- !is.na(observed_LODs)
      p_values <- rep(NA_real_, length(observed_LODs))
      num_perm_lods <- length(perm_LODs_all)

      if (num_perm_lods > 0) {
        p_values[valid_obs_indices] <- sapply(observed_LODs[valid_obs_indices], function(obs_lod) {
          count_ge <- sum(perm_LODs_all >= obs_lod, na.rm = TRUE)
          (count_ge + 1) / (num_perm_lods + 1)
        })
      } else {
        warning("Permutation resulted in zero valid LOD scores. Cannot calculate p-values.", call.=FALSE)
      }

      q_values <- rep(NA_real_, length(p_values))
      valid_p_indices <- !is.na(p_values)
      if(any(valid_p_indices)) {
        q_values[valid_p_indices] <- stats::p.adjust(p_values[valid_p_indices], method = "BH")
      }

      results_df_genome$p_value <- p_values
      results_df_genome$q_value <- q_values

      # Calculate the LOD threshold corresponding to the q-value cutoff
      significant_lods <- results_df_genome$LOD[!is.na(results_df_genome$q_value) & results_df_genome$q_value <= q_threshold]
      if(length(significant_lods) > 0) {
        fdr_lod_threshold <- min(significant_lods, na.rm = TRUE)
        msg(paste0("LOD threshold corresponding to q-value <= ", q_threshold, ": ", round(fdr_lod_threshold, 3)))
      } else {
        msg(paste0("No bins found with q-value <= ", q_threshold, "."))
        fdr_lod_threshold <- Inf # Set to Inf if none are significant
      }

    } else {
      msg("Permutation testing did not yield valid results. Skipping p-value/q-value calculation.")
      results_df_genome$p_value <- NA_real_
      results_df_genome$q_value <- NA_real_
    }

  } else {
    msg("Skipping permutation test (nperm = 0). p-values and q-values will not be calculated.")
    results_df_genome$p_value <- NA_real_
    results_df_genome$q_value <- NA_real_
  }

  # --- Gene Annotation ---
  genes_annotated <- NULL
  # Perform annotation if DBs are available
  if (!is.null(loaded_txdb) && !is.null(loaded_orgdb)) {
    bins_for_annotation <- NULL
    if (nperm > 0 && "q_value" %in% names(results_df_genome)) {
      significant_bins <- !is.na(results_df_genome$q_value) & results_df_genome$q_value <= q_threshold
      if (any(significant_bins)) {
        msg("Annotating genes in regions with q-value <= ", q_threshold, "...")
        bins_for_annotation <- results_df_genome[significant_bins, ]
      } else {
        msg("No bins found below q-value threshold ", q_threshold, ". No genes annotated based on significance.")
      }
    } else if (nperm == 0) {
      msg("Annotating all genes overlapping analyzed bins (nperm=0)...")
      bins_for_annotation <- results_df_genome
    }

    if (!is.null(bins_for_annotation) && nrow(bins_for_annotation) > 0) {
      genes_annotated <- internal_annotate_qtl(
        results_df = bins_for_annotation,
        q_threshold = q_threshold, # Pass for context
        txdb = loaded_txdb,
        orgdb = loaded_orgdb,
        common_chroms = chrom_order
      )
      if (!is.null(genes_annotated)) {
        msg(paste("Found", nrow(genes_annotated), "unique genes overlapping selected regions."))
      }
    } else if (nperm > 0) {
      # Message already printed if no bins were significant
    }

  } else {
    msg("Skipping annotation: txdb or orgdb object not available.")
  }


  # --- Plotting ---
  p_genome <- NULL
  if (plot) {
    msg("Generating genome-wide plot...")
    if (nrow(results_df_genome) > 0) {
      p_genome <- internal_plot_genome_wide(
        results_df = results_df_genome,
        fdr_lod_threshold = fdr_lod_threshold, # Pass the calculated LOD threshold
        genes_annotated = genes_annotated,
        nperm = nperm,
        q_threshold = q_threshold,
        mode = mode,
        chrom_order = chrom_order
      )
      if (interactive() && !is.null(p_genome)) print(p_genome)
    } else {
      msg("Skipping plot generation: No results to plot.")
    }
  }

  # --- Return Results ---
  msg("Analysis complete.")
  invisible(list(
    results = results_df_genome,
    fdr_lod_threshold = fdr_lod_threshold, # Return the threshold value
    plot = p_genome,
    annotated_genes = genes_annotated,
    parameters = list(
      N = N, mode = mode, res = res, cM = cM, filter = filter,
      nperm = nperm, q_threshold = q_threshold,
      seed = seed, parallel = parallel,
      p_recombination = p_rec
    )
  ))
}


# ---- Internal Data Loading ----
# (No changes needed in internal_load_two_pools, internal_bin_chromosome_data, internal_unify_bins)
#' Load and preprocess data for two pools
#' @keywords internal
internal_load_two_pools <- function(pool1, pool2, res, filter,
                                    col_chr ,
                                    col_pos  ,
                                    col_a    ,
                                    col_b  ,
                                    assume_chr ,
                                    header) {

  read_pool_data <- function(pool_input,
                             col_chr=col_chr,
                             col_pos  ,
                             col_a    ,
                             col_b  ,
                             assume_chr ,
                             header) {

    # --- load -------------------------------------------------------
    if (is.character(pool_input)) {
      df <- read.table(pool_input, header = header, stringsAsFactors = FALSE)
    } else if (is.data.frame(pool_input)) {
      df <- pool_input
      if (!header) names(df) <- NULL
    } else {
      stop("Input pool must be a file path or data frame.")
    }

    # --- assign column names when header = FALSE -------------------
    if (!header) {
      n <- ncol(df)
      if (n == 3)       names(df) <- c(col_pos, col_a, col_b)
      else if (n == 4)  names(df) <- c(col_chr, col_pos, col_a, col_b)
      else stop("File with header = FALSE must have 3 or 4 columns.")
    }

    # --- rename user-supplied column names to internal a/b ----------
    req <- c(col_pos, col_a, col_b)
    if (is.null(col_chr) || !(col_chr %in% names(df))) {
      df[[col_chr]] <- assume_chr         # single chromosome fallback
    }

    if (!is.null(col_chr)) req <- c(col_chr, req)
    if (!all(req %in% names(df)))
      stop("Missing required columns: ", paste(setdiff(req, names(df)), collapse = ", "))

    names(df)[match(c(col_chr, col_pos, col_a, col_b), names(df))] <-
      c("chr", "pos", "a", "b")        # canonical names

    df
  }

  df1 <- read_pool_data(pool1, col_chr = col_chr, col_pos = col_pos,
                        col_a = col_a, col_b = col_b, assume_chr = assume_chr,
                        header = header)
  df2 <- read_pool_data(pool2, col_chr = col_chr, col_pos = col_pos,
                        col_a = col_a, col_b = col_b, assume_chr = assume_chr,
                        header = header)
  if (nrow(df1) == 0 || nrow(df2) == 0) {
    stop("One or both input pools have no valid data rows after initial checks.")
  }
  common_chroms <- sort(intersect(unique(df1$chr), unique(df2$chr)))
  if (length(common_chroms) == 0) {
    stop("No common chromosomes found between the two input pools.")
  }
  msg(paste("Found", length(common_chroms), "common chromosomes:", paste(common_chroms, collapse=", ")))
  chrom_data_list <- list()
  orig_bins_list <- list()
  for (chrom in common_chroms) {
    df1_chr <- df1[df1$chr == chrom, ]
    df2_chr <- df2[df2$chr == chrom, ]
    if (nrow(df1_chr) == 0 || nrow(df2_chr) == 0) {
      warning("Skipping chromosome ", chrom, " as it has no data in one or both pools after filtering by chromosome.", call.=FALSE)
      next
    }
    binned1 <- internal_bin_chromosome_data(df1_chr, res, filter, pool_name = "Pool 1", chrom_name = chrom)
    binned2 <- internal_bin_chromosome_data(df2_chr, res, filter, pool_name = "Pool 2", chrom_name = chrom)
    if (is.null(binned1) || is.null(binned2)) {
      warning("Skipping chromosome ", chrom, " due to issues during binning (e.g., no data left after filtering).", call.=FALSE)
      next
    }
    min_start <- min(binned1$bins[1], binned2$bins[1])
    max_end <- max(tail(binned1$bins, 1), tail(binned2$bins, 1))
    full_bins <- seq(min_start, max_end, by = res)
    if (length(full_bins) < 2) {
      warning("Skipping chromosome ", chrom, ": cannot create valid bins (min_start=", min_start, ", max_end=", max_end, ", res=", res, ").", call.=FALSE)
      next
    }
    unify1 <- internal_unify_bins(binned1$means, binned1$variances, binned1$counts, binned1$bins, full_bins)
    unify2 <- internal_unify_bins(binned2$means, binned2$variances, binned2$counts, binned2$bins, full_bins)
    T_chr <- length(unify1$means)
    if (T_chr == 0) {
      warning("Skipping chromosome ", chrom, " as it has no bins after unifying.", call.=FALSE)
      next
    }
    chrom_data_list[[chrom]] <- list(
      y1 = unify1$means, y_var1 = unify1$variances, d1 = unify1$counts,
      y2 = unify2$means, y_var2 = unify2$variances, d2 = unify2$counts,
      T_ = T_chr,
      bins = unify1$bins
    )
    orig_bins_list[[chrom]] <- list(pool1 = binned1$bins, pool2 = binned2$bins)
    msg(paste("  Chromosome", chrom, ":",
                  sum(unify1$counts), "reads (Pool 1),",
                  sum(unify2$counts), "reads (Pool 2) across",
                  T_chr, "bins."))
  }
  return(list(chrom_data = chrom_data_list,
              orig_bins = orig_bins_list,
              common_chroms = common_chroms))
}

#' Bin data for a single chromosome from one pool
#' @keywords internal
internal_bin_chromosome_data <- function(df_chr, binsize, filt = TRUE, pool_name = "", chrom_name = "") {
  df_chr$pos <- as.numeric(df_chr$pos)
  df_chr$a <- as.numeric(df_chr$a)
  df_chr$b <- as.numeric(df_chr$b)
  df_chr <- df_chr[!is.na(df_chr$pos) & !is.na(df_chr$a) & !is.na(df_chr$b), ]
  if (nrow(df_chr) == 0) {
    msg("  No valid markers for ", pool_name, " on chromosome ", chrom_name, " before binning.")
    return(NULL)
  }
  n_before_filt_fix <- nrow(df_chr)
  if (filt) {
    df_chr <- df_chr[df_chr$a > 0 & df_chr$b > 0, ]
    n_after_filt_fix <- nrow(df_chr)
    if (n_before_filt_fix > n_after_filt_fix) {
      msg(paste("  Filtered", n_before_filt_fix - n_after_filt_fix,
                    "fixated markers (a=0 or b=0) for", pool_name, "on chr", chrom_name))
    }
  }
  if (nrow(df_chr) == 0) {
    msg("  No markers left for ", pool_name, " on chromosome ", chrom_name, " after fixation filter.")
    return(NULL)
  }
  bin_starts <- floor(df_chr$pos / binsize) * binsize
  min_pos <- min(df_chr$pos)
  max_pos <- max(df_chr$pos)
  global_bin_edges <- seq(floor(min_pos / binsize) * binsize,
                          ceiling(max_pos / binsize) * binsize,
                          by = binsize)
  if (length(global_bin_edges) < 2) {
    global_bin_edges <- c(floor(min_pos / binsize) * binsize,
                          floor(min_pos / binsize) * binsize + binsize)
  }
  agg_a <- stats::aggregate(df_chr$a, by = list(bin_start = bin_starts), FUN = sum)
  agg_b <- stats::aggregate(df_chr$b, by = list(bin_start = bin_starts), FUN = sum)
  binned_counts <- merge(agg_a, agg_b, by = "bin_start", all = TRUE)
  names(binned_counts) <- c("bin_start", "a_sum", "b_sum")
  binned_counts$a_sum[is.na(binned_counts$a_sum)] <- 0
  binned_counts$b_sum[is.na(binned_counts$b_sum)] <- 0
  binned_counts$total_count <- binned_counts$a_sum + binned_counts$b_sum
  n_bins_before_filt_cov <- nrow(binned_counts)
  if (filt && nrow(binned_counts) > 1) {
    total_counts <- binned_counts$total_count
    med_ <- stats::median(total_counts)
    mad_ <- stats::median(abs(total_counts - med_))
    if (mad_ == 0 && med_ > 0) {
      mad_ <- med_ * 0.1
      cutoff <- med_ + 20 * mad_
      msg(paste("  Coverage MAD is 0 for", pool_name, "on chr", chrom_name,
                    ". Using cutoff based on 10% of median coverage:", round(cutoff,1)))
    } else if (mad_ > 0) {
      cutoff <- med_ + 20 * mad_
      msg(paste("  Coverage filter cutoff (median + 20*MAD) for", pool_name,
                    "on chr", chrom_name, ":", round(cutoff, 1),
                    "(median=", round(med_,1), ", MAD=", round(mad_,1), ")"))
    } else {
      cutoff <- Inf
    }
    binned_counts_filtered <- binned_counts[total_counts <= cutoff, ]
    n_bins_after_filt_cov <- nrow(binned_counts_filtered)
    if (n_bins_before_filt_cov > n_bins_after_filt_cov) {
      msg(paste("  Filtered", n_bins_before_filt_cov - n_bins_after_filt_cov,
                    "bins with coverage >", round(cutoff, 1), "for", pool_name, "on chr", chrom_name))
    }
    binned_counts <- binned_counts_filtered
  }
  if (nrow(binned_counts) == 0) {
    msg("  No bins left for ", pool_name, " on chromosome ", chrom_name, " after coverage filter.")
    return(NULL)
  }
  num_bins <- length(global_bin_edges) - 1
  means_out <- numeric(num_bins)
  variances_out <- rep(Inf, num_bins)
  counts_out <- numeric(num_bins)
  bin_indices <- match(binned_counts$bin_start, global_bin_edges[-length(global_bin_edges)])
  good <- which(!is.na(bin_indices) & bin_indices >= 1 & bin_indices <= num_bins)
  if (length(good) == 0) {
    msg("  No bins map cleanly onto global grid for ", pool_name,
            " on chromosome ", chrom_name, ".")
    return(NULL)
  }
  bin_indices   <- bin_indices[good]
  binned_counts <- binned_counts[good, , drop = FALSE]
  valid_indices <- !is.na(bin_indices)
  target_indices <- bin_indices[valid_indices]
  if (any(target_indices > num_bins | target_indices < 1)) {
    warning("Bin index calculation error on chr ", chrom_name, ". Some indices are out of expected range [1, ", num_bins, "].", call. = FALSE)
    invalid_target <- target_indices > num_bins | target_indices < 1
    target_indices <- target_indices[!invalid_target]
    binned_counts <- binned_counts[valid_indices, ][!invalid_target, ]
    if(length(target_indices) == 0) return(NULL)
  }
  means_out[target_indices] <- binned_counts$a_sum
  counts_out[target_indices] <- binned_counts$total_count
  valid_counts_idx_in_binned <- which(binned_counts$total_count > 0)
  target_indices_valid_counts <- target_indices[valid_counts_idx_in_binned]
  if (length(target_indices_valid_counts) > 0) {
    a_sum_valid <- binned_counts$a_sum[valid_counts_idx_in_binned]
    total_count_valid <- binned_counts$total_count[valid_counts_idx_in_binned]
    p_hat <- a_sum_valid / total_count_valid
    p_hat <- pmax(1e-9, pmin(1 - 1e-9, p_hat))
    calc_variances <- total_count_valid * p_hat * (1 - p_hat)
    variances_out[target_indices_valid_counts] <- pmax(1e-9, calc_variances)
  }
  list(
    means = means_out,
    variances = variances_out,
    counts = counts_out,
    bins = global_bin_edges
  )
}

#' Unify bin data to a common set of bin edges
#' @keywords internal
internal_unify_bins <- function(means, variances, counts, old_bins, new_bins) {
  nb <- length(new_bins) - 1
  out_means <- numeric(nb)
  out_vars <- rep(Inf, nb)
  out_counts <- numeric(nb)
  old_starts <- old_bins[-length(old_bins)]
  new_starts <- new_bins[-length(new_bins)]
  match_indices <- match(old_starts, new_starts, nomatch = 0)
  valid_old_bin_indices <- which(match_indices != 0)
  target_new_bin_indices <- match_indices[valid_old_bin_indices]
  out_means[target_new_bin_indices] <- means[valid_old_bin_indices]
  out_vars[target_new_bin_indices] <- variances[valid_old_bin_indices]
  out_counts[target_new_bin_indices] <- counts[valid_old_bin_indices]
  list(
    means = out_means,
    variances = out_vars,
    counts = out_counts,
    bins = new_bins
  )
}


# ---- Internal Core Computation ----

#' Perform core multipool computation using Kalman filter
#' @keywords internal
internal_do_computation <- function(y1, y_var1, d1,
                                    y2, y_var2, d2,
                                    T_, N, p,
                                    mode) {
  k1 <- internal_kalman(y = y1, y_var = y_var1, d = d1, T_ = T_, N = N, p = p)
  k2 <- internal_kalman(y = y2, y_var = y_var2, d = d2, T_ = T_, N = N, p = p)
  lod_calc1 <- internal_calc_lods_multicoupled(
    mu_pstr_list = list(k1$mu_pstr),
    V_pstr_list = list(k1$V_pstr),
    T_ = T_, N = N
  )
  lod_calc2 <- internal_calc_lods_multicoupled(
    mu_pstr_list = list(k2$mu_pstr),
    V_pstr_list = list(k2$V_pstr),
    T_ = T_, N = N
  )
  lod_calc_joint <- internal_calc_lods_multicoupled(
    mu_pstr_list = list(k1$mu_pstr, k2$mu_pstr),
    V_pstr_list = list(k1$V_pstr, k2$V_pstr),
    T_ = T_, N = N
  )
  if (mode == "replicates") {
    LOD <- lod_calc_joint$LOD
    mu_MLE <- lod_calc_joint$mu_MLE
  } else {
    LOD <- lod_calc1$LOD + lod_calc2$LOD - lod_calc_joint$LOD
    LOD[LOD < 0] <- 0
    mu_MLE <- lod_calc1$mu_MLE - lod_calc2$mu_MLE
  }
  list(
    LOD = LOD,
    mu_MLE = mu_MLE,
    mu_pstr1 = k1$mu_pstr,
    V_pstr1 = k1$V_pstr,
    mu_pstr2 = k2$mu_pstr,
    V_pstr2 = k2$V_pstr
  )
}

#' Kalman Filter and Smoother for Allele Counts
#' @keywords internal
internal_kalman <- function(y, y_var, d, T_, N, p) {
  mu_pred <- numeric(T_)
  V_pred <- numeric(T_)
  mu_filt <- numeric(T_)
  V_filt <- numeric(T_)
  K <- numeric(T_)
  c_lik <- numeric(T_)
  mu_pstr <- numeric(T_)
  V_pstr <- numeric(T_)
  J <- numeric(T_)
  mu_init <- 0.5 * N
  V_init <- 0.25 * N
  A <- 1 - 2 * p
  C <- d / N
  S <- p * (1 - p) * N
  mu_pred[1] <- mu_init
  V_pred[1] <- max(V_init, 1e-9)
  if (is.infinite(y_var[1]) || C[1] == 0) {
    K[1] <- 0
    c_lik[1] <- 1
    mu_filt[1] <- mu_pred[1]
    V_filt[1] <- V_pred[1]
  } else {
    innov_var <- C[1]^2 * V_pred[1] + y_var[1]
    if (innov_var <= 1e-9) {
      warning("Kalman Filter: Non-positive innovation variance at t=1. Setting K=0.", call.=FALSE)
      K[1] <- 0
      c_lik[1] <- 1
      mu_filt[1] <- mu_pred[1]
      V_filt[1] <- V_pred[1]
    } else {
      K[1] <- V_pred[1] * C[1] / innov_var
      innov_mean <- C[1] * mu_pred[1]
      c_lik[1] <- stats::dnorm(y[1], mean = innov_mean, sd = sqrt(innov_var))
      c_lik[1] <- max(c_lik[1], 1e-300)
      mu_filt[1] <- mu_pred[1] + K[1] * (y[1] - C[1] * mu_pred[1])
      V_filt[1] <- (1 - K[1] * C[1]) * V_pred[1]
    }
  }
  V_filt[1] <- max(V_filt[1], 1e-9)
  for (t in 2:T_) {
    mu_pred[t] <- A * mu_filt[t - 1] + p * N
    V_pred[t] <- A^2 * V_filt[t - 1] + S
    V_pred[t] <- max(V_pred[t], 1e-9)
    if (is.infinite(y_var[t]) || C[t] == 0) {
      K[t] <- 0
      c_lik[t] <- 1
      mu_filt[t] <- mu_pred[t]
      V_filt[t] <- V_pred[t]
    } else {
      innov_var <- C[t]^2 * V_pred[t] + y_var[t]
      if (innov_var <= 1e-9) {
        warning(paste("Kalman Filter: Non-positive innovation variance at t=", t, ". Setting K=0."), call.=FALSE)
        K[t] <- 0
        c_lik[t] <- 1
        mu_filt[t] <- mu_pred[t]
        V_filt[t] <- V_pred[t]
      } else {
        K[t] <- V_pred[t] * C[t] / innov_var
        innov_mean <- C[t] * mu_pred[t]
        c_lik[t] <- stats::dnorm(y[t], mean = innov_mean, sd = sqrt(innov_var))
        c_lik[t] <- max(c_lik[t], 1e-300)
        mu_filt[t] <- mu_pred[t] + K[t] * (y[t] - C[t] * mu_pred[t])
        V_filt[t] <- (1 - K[t] * C[t]) * V_pred[t]
      }
    }
    V_filt[t] <- max(V_filt[t], 1e-9)
  }
  mu_pstr[T_] <- mu_filt[T_]
  V_pstr[T_] <- V_filt[T_]
  for (t in seq(T_ - 1, 1, by = -1)) {
    J[t] <- V_filt[t] * A / V_pred[t + 1]
    mu_pstr[t] <- mu_filt[t] + J[t] * (mu_pstr[t + 1] - mu_pred[t + 1])
    V_pstr[t] <- V_filt[t] + J[t]^2 * (V_pstr[t + 1] - V_pred[t + 1])
    V_pstr[t] <- max(V_pstr[t], 1e-9)
  }
  logLik <- sum(log(c_lik))
  list(
    mu_pstr = mu_pstr,
    V_pstr = V_pstr,
    logLik = logLik
  )
}

# Calculate LOD Scores using Optimization (R Implementation)
# @keywords internal

#' Calculate LOD Scores using Optimization (R Implementation - Python Logic, Optimized)
#'
#' Calculates LOD scores based on the posterior estimates from the Kalman filter.
#' This version replicates the specific calculation logic found in the original
#' Python implementation's `calcLODs_multicoupled` function but uses matrix
#' operations for significantly improved performance compared to the direct loop-based
#' translation.
#'
#' @param mu_pstr_list List of posterior mean vectors (one per pool).
#' @param V_pstr_list List of posterior variance vectors (one per pool).
#' @param T_ Integer, number of bins.
#' @param N Numeric, effective population size.
#'
#' @return A list containing `LOD` (vector of LOD scores) and `mu_MLE` (vector
#'   of maximum likelihood estimates of the underlying count N*p_alt).
#' @keywords internal
internal_calc_lods_multicoupled <- function(mu_pstr_list,
                                            V_pstr_list,
                                            T_,
                                            N) {

  # --- Input Validation ---
  if (length(mu_pstr_list) != length(V_pstr_list) ||
      any(vapply(mu_pstr_list, length, 0L) != T_))
    stop("Posterior mean/var lists must have equal length and length T_")

  num_pools <- length(mu_pstr_list)
  if (num_pools == 0) stop("Empty posterior lists provided.")

  # --- Setup Grid and Constants (matching Python) ---
  delta <- DELTA_X # Use the globally defined DELTA_X
  p_alt_grid <- seq(delta, 1.0 - delta, by = delta) # Grid for p' (allele frequency)
  cnt_grid <- N * p_alt_grid                        # Grid for counts j = N*p'
  grid_size <- length(p_alt_grid)

  # --- Precompute P(x=j | p=p') matrix (p_precomp) ---
  # Rows: alternative p' (p_alt_idx), Columns: count grid j (cnt_grid_idx)
  # Using outer to vectorize the calculation of the matrix
  p_alt_means <- N * p_alt_grid
  p_alt_vars <- N * p_alt_grid * (1.0 - p_alt_grid)
  p_alt_sds <- sqrt(pmax(p_alt_vars, 1e-9))

  # Calculate the normal PDF for all combinations of cnt_grid and p_alt parameters
  # outer will apply dnorm to each element of cnt_grid vs each p_alt mean/sd
  # Need to transpose because outer iterates the second argument faster
  p_precomp <- t(outer(cnt_grid, 1:grid_size, FUN = function(cnt, idx) {
    stats::dnorm(cnt, mean = p_alt_means[idx], sd = p_alt_sds[idx])
  }))
  # Resulting p_precomp matrix has dimensions: grid_size x grid_size
  # Row i corresponds to p_alt_grid[i], Column j corresponds to cnt_grid[j]

  # --- Precompute log P(x=j) under null (logreweighter) ---
  mu_initial <- 0.5 * N
  V_initial <- 0.25 * N
  sd_initial <- sqrt(max(V_initial, 1e-9))
  logreweighter <- internal_lognormpdf(cnt_grid, mu_initial, sd_initial) # Vector length grid_size

  # --- Initialize Output Vectors ---
  LOD <- numeric(T_)
  mu_MLE <- numeric(T_) # Stores N * p_alt_MLE

  # --- Loop through each bin (time point t) ---
  # This outer loop remains as posteriors change per bin
  for (t in 1:T_) {
    # Extract posterior means and variances for *all pools* at the current bin t
    mu_pstr_t_vec <- vapply(mu_pstr_list, `[`, numeric(1), t)
    V_pstr_t_vec <- vapply(V_pstr_list, `[`, numeric(1), t)
    sd_pstr_t_vec <- sqrt(pmax(V_pstr_t_vec, 1e-9))

    # --- Vectorized calculation for bin t ---

    # 1. Calculate logtemp matrix (grid_size x num_pools)
    #    logtemp_kj = log( P(x=j | y_k,t) / P(x=j | null) )
    #    Use outer to apply internal_lognormpdf across cnt_grid and pool parameters
    logtemp_matrix <- outer(cnt_grid, 1:num_pools, FUN = function(cnt, k) {
      internal_lognormpdf(cnt, mu_pstr_t_vec[k], sd_pstr_t_vec[k])
    }) - logreweighter # Recycle logreweighter vector across columns

    # 2. Calculate scalers (vector of length num_pools)
    #    Find max logtemp for each pool (column) to prevent underflow/overflow
    scalers <- apply(logtemp_matrix, 2, max, na.rm = TRUE)
    scalers[!is.finite(scalers)] <- -700.0 # Handle cases where all are -Inf

    # 3. Calculate scaled exponential matrix (grid_size x num_pools)
    #    exp(logtemp_kj - scaler_k)
    #    Expand scalers to match matrix dimensions for subtraction
    exp_logtemp_scaled_matrix <- exp(logtemp_matrix - matrix(scalers, nrow = grid_size, ncol = num_pools, byrow = TRUE))

    # 4. Calculate dot products matrix (grid_size x num_pools)
    #    dotprod_ik = sum_j [ P(x=j | p'_i) * exp(logtemp_kj - scaler_k) ]
    #    This is matrix multiplication: p_precomp %*% exp_logtemp_scaled_matrix
    dot_prod_matrix <- p_precomp %*% exp_logtemp_scaled_matrix
    # Resulting matrix has dimensions: grid_size x num_pools

    # 5. Calculate log of dot products + scaler (grid_size x num_pools)
    #    log(max(dotprod_ik, 1e-300)) + scaler_k
    log_dot_prod_matrix <- log(pmax(dot_prod_matrix, 1e-300)) + matrix(scalers, nrow = grid_size, ncol = num_pools, byrow = TRUE)

    # 6. Sum logs across pools (vector of length grid_size)
    #    logallsums_i = sum_k [ log(dotprod_ik) + scaler_k ]
    logallsums_grid_t <- rowSums(log_dot_prod_matrix, na.rm = TRUE) # Sum across pools (columns)

    # --- Find Maximum and Calculate LOD for bin t ---
    max_log_sum_t <- max(logallsums_grid_t, na.rm = TRUE)
    max_idx_t <- which.max(logallsums_grid_t)
    if(length(max_idx_t) > 1) max_idx_t <- max_idx_t[1] # Handle ties

    mu_MLE[t] <- N * p_alt_grid[max_idx_t]
    LOD[t] <- log10(N) + log10(delta) + max_log_sum_t / log(10.0)

  } # End loop over bins t

  # Final cleanup
  LOD[LOD < 0 | !is.finite(LOD)] <- 0

  list(LOD = LOD, mu_MLE = mu_MLE)
}

# Make sure internal_lognormpdf is defined correctly as before:
#' Safe Log Normal PDF Calculation (R version)
#' @keywords internal
internal_lognormpdf <- function(x, mu, sigma) {
  # Ensure sigma is positive
  sigma <- pmax(sigma, 1e-9)
  # Use dnorm with log = TRUE
  log_pdf <- stats::dnorm(x, mean = mu, sd = sigma, log = TRUE)
  # Clamp extreme negative values to avoid -Inf issues downstream
  log_pdf[is.infinite(log_pdf) & log_pdf < 0] <- -700 # Approx log(1e-300)
  log_pdf[!is.finite(log_pdf)] <- -700 # Handle NaN or +Inf just in case
  return(log_pdf)
}


# internal_calc_lods_multicoupled <- function(mu_pstr_list,
#                                             V_pstr_list,
#                                             T_,
#                                             N)
# {
#   if (length(mu_pstr_list) != length(V_pstr_list) ||
#       any(vapply(mu_pstr_list, length, 0L) != T_))
#     stop("Posterior mean/var lists must have equal length and length T_")
#
#   G     <- seq(DELTA_X, 1 - DELTA_X, by = DELTA_X)        # freq grid
#   cnt   <- G * N
#   Gsize <- length(G)
#
#   # constant importance density (log space)
#   base_lp <- stats::dnorm(cnt,
#                           mean = 0.5 * N,
#                           sd   = sqrt(0.25 * N),
#                           log  = TRUE)
#
#   # accumulator grid: rows = grid, cols = bins
#   lp_sum <- matrix(0, nrow = Gsize, ncol = T_)
#
#   for (k in seq_along(mu_pstr_list)) {
#     mu  <- mu_pstr_list[[k]]
#     sig <- sqrt(pmax(V_pstr_list[[k]], 1e-9))
#
#     # Build mu & sigma matrices by row repetition (fast in C)
#     mu_mat  <- matrix(rep(mu, each = Gsize),  nrow = Gsize)
#     sig_mat <- matrix(rep(sig, each = Gsize), nrow = Gsize)
#
#     lp_k <- stats::dnorm(cnt, mu_mat, sig_mat, log = TRUE) - base_lp
#     lp_sum <- lp_sum + lp_k
#
#   }
#
#   # For each column (bin), find maximum log-likelihood ratio and its index
#   max_log <- matrixStats::colMaxs(lp_sum)
#   max_idx <- max.col(t(lp_sum), ties.method = "first")           # gives first max pos
#
#   p_alt   <- G[max_idx]
#   LOD     <- log10(N) + log10(DELTA_X) + max_log / log(10)
#   LOD[LOD < 0 | !is.finite(LOD)] <- 0
#
#   LOD <- pmin(LOD, 300)               # hard ceiling as in legacy code
#
#   list(LOD = LOD,
#        mu_MLE = p_alt * N)
# }


#' Safe Log Normal PDF Calculation (R version)
#' @keywords internal
internal_lognormpdf <- function(x, mu, sigma) {
  log_pdf <- numeric(length(x))
  sigma_threshold <- 1e-9
  valid_sigma <- sigma > sigma_threshold
  if (any(valid_sigma)) {
    sigma_valid <- if(length(sigma) == 1) rep(sigma, sum(valid_sigma)) else sigma[valid_sigma]
    mu_valid <- if(length(mu) == 1) rep(mu, sum(valid_sigma)) else mu[valid_sigma]
    x_valid <- x[valid_sigma]
    log_pdf[valid_sigma] <- stats::dnorm(x_valid, mean = mu_valid, sd = sigma_valid, log = TRUE)
  }
  if (any(!valid_sigma)) {
    is_close <- abs(x[!valid_sigma] - mu) < 1e-6
    log_pdf[!valid_sigma & is_close] <- log(1e300)
    log_pdf[!valid_sigma & !is_close] <- log(1e-300)
  }
  log_pdf[is.infinite(log_pdf) & log_pdf < 0] <- log(1e-300)
  log_pdf[!is.finite(log_pdf)] <- log(1e-300)
  return(log_pdf)
}


# ---- Internal Permutation Test ----

#' Permutation test for Multipool (Internal Implementation - FDR/q-value)
#'
#' Performs permutation testing by shuffling the *total read counts* within each bin
#' between the two pools across the entire genome. This disrupts the association
#' between the pools at each genomic location while preserving the overall
#' distribution of read counts.
#'
#' @param chrom_data_list List of chromosome data (output from internal_load_two_pools$chrom_data).
#' @param N Effective number of individuals.
#' @param p Recombination rate.
#' @param mode "replicates" or "contrast".
#' @param nperm Number of permutations.
#' @return A numeric vector containing all LOD scores from all bins across all permutations.
#' @keywords internal
internal_multipool_permutation_fdr <- function(chrom_data_list, N, p, mode, nperm,
                                               seed     = NULL,
                                               parallel = FALSE) {

  if (!is.null(seed)) set.seed(seed)

  # Filter out chromosomes with T_ < 2 before proceeding
  valid_chrom_data <- Filter(function(cd) !is.null(cd) && cd$T_ >= 2, chrom_data_list)
  if (length(valid_chrom_data) == 0) {
    warning("Permutation test skipped: No chromosomes with sufficient bins (>= 2) found.", call. = FALSE)
    return(numeric(0)) # Return empty numeric vector
  }

  # Calculate total number of bins across valid chromosomes and store original data
  chrom_lengths <- sapply(valid_chrom_data, `[[`, "T_")
  total_bins <- sum(chrom_lengths)
  if (total_bins == 0) {
    warning("Permutation test skipped: No bins found across valid chromosomes.", call. = FALSE)
    return(numeric(0))
  }

  original_y1 <- unlist(lapply(valid_chrom_data, `[[`, "y1"))
  original_y_var1 <- unlist(lapply(valid_chrom_data, `[[`, "y_var1"))
  original_d1 <- unlist(lapply(valid_chrom_data, `[[`, "d1"))
  original_y2 <- unlist(lapply(valid_chrom_data, `[[`, "y2"))
  original_y_var2 <- unlist(lapply(valid_chrom_data, `[[`, "y_var2"))
  original_d2 <- unlist(lapply(valid_chrom_data, `[[`, "d2"))

  # Pre-allocate list to store LOD vectors from each permutation
  perm_LODs_list <- vector("list", nperm)

  # --- Permutation Loop ---
  # for (i in seq_len(nperm)) {
    # if (i %% 50 == 0 || i == nperm) msg(paste("  Permutation", i, "of", nperm))
  run_single_perm <- function(i) {
    if (!is.null(seed)) set.seed(seed + i)

    # Create permuted data vectors
    permuted_y1 <- original_y1
    permuted_y_var1 <- original_y_var1
    permuted_d1 <- original_d1
    permuted_y2 <- original_y2
    permuted_y_var2 <- original_y_var2
    permuted_d2 <- original_d2

    # Shuffle counts between pools within each bin
    for (j in seq_len(total_bins)) {
      if (runif(1) < 0.5) { # 50% chance of swapping
        # Swap y (allele A counts)
        temp_y1 <- permuted_y1[j]
        permuted_y1[j] <- permuted_y2[j]
        permuted_y2[j] <- temp_y1
        # Swap d (total counts)
        temp_d1 <- permuted_d1[j]
        permuted_d1[j] <- permuted_d2[j]
        permuted_d2[j] <- temp_d1
        # Note: y_var is a function of y and d, so we should recalculate it.
        # However, for simplicity and consistency with the original structure,
        # we will keep the original y_var values. A more rigorous approach might involve
        # recalculating y_var based on the permuted y and d.
      }
    }

    # --- Re-run computation on shuffled data (per chromosome) ---
    lod_scores_this_perm <- numeric(total_bins) # Vector to store all LODs for this perm
    current_pos <- 0
    for (chrom_idx in seq_along(chrom_lengths)) {
      T_chr <- chrom_lengths[chrom_idx]
      chr_indices <- (current_pos + 1):(current_pos + T_chr)

      # Extract shuffled data for this chromosome
      y1p_chr <- permuted_y1[chr_indices]
      yv1p_chr <- permuted_y_var1[chr_indices]
      d1p_chr <- permuted_d1[chr_indices]
      y2p_chr <- permuted_y2[chr_indices]
      yv2p_chr <- permuted_y_var2[chr_indices]
      d2p_chr <- permuted_d2[chr_indices]

      # Run computation
      comp_perm <- internal_do_computation(
        y1 = y1p_chr, y_var1 = yv1p_chr, d1 = d1p_chr,
        y2 = y2p_chr, y_var2 = yv2p_chr, d2 = d2p_chr,
        T_ = T_chr, N = N, p = p, mode = mode
      )

      # Store all LOD scores for this chromosome
      lod_scores_this_perm[chr_indices] <- comp_perm$LOD
      current_pos <- current_pos + T_chr
    } # End chromosome loop

    # Store the vector of all LOD scores for this permutation
    perm_LODs_list[[i]] <- lod_scores_this_perm

  } # End permutation loop

  perm_LODs_list <- if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    future.apply::future_lapply(seq_len(nperm), run_single_perm,
                                future.seed = TRUE)
  } else {
    if (parallel)
      msg("future.apply not available – running permutations serially.")
    lapply(seq_len(nperm), run_single_perm)
  }

  # Combine all permuted LOD scores into one large vector
  perm_LODs_all <- unlist(perm_LODs_list)
  # Remove any NAs or non-finite values that might have occurred
  perm_LODs_all <- perm_LODs_all[is.finite(perm_LODs_all)]

  if (length(perm_LODs_all) == 0) {
    warning("No finite LOD scores generated during permutations.", call. = FALSE)
  }

  return(perm_LODs_all)
}


# ---- Internal Annotation of Yeast Genes ----

#' Annotate QTL regions with yeast genes (Internal Implementation)
#'
#' Finds genes overlapping bins. If a q_threshold is provided (implicitly via
#' results_df containing a q_value column), only genes overlapping bins
#' meeting the threshold are returned. Otherwise (if nperm=0), all
#' overlapping genes are returned.
#'
#' @param results_df Data frame of bins to annotate (pre-filtered based on significance or all bins if nperm=0).
#'           Must include chr, pos, LOD.
#' @param q_threshold Optional numeric q-value significance threshold (only used for messaging).
#' @param txdb TxDb object (e.g., TxDb.Scerevisiae.UCSC.sacCer3.sgdGene).
#' @param orgdb OrgDb object (e.g., org.Sc.sgd.db).
#' @param common_chroms Character vector of common chromosome names used in analysis (for ordering).
#' @return Data frame of unique annotated genes, or NULL.
#' @keywords internal
internal_annotate_qtl <- function(results_df, q_threshold = NULL, txdb, orgdb, common_chroms) { # Corrected argument name

  # Input checks
  if (is.null(txdb) || is.null(orgdb)) {
    msg("Annotation skipped: Missing txdb or orgdb.")
    return(NULL)
  }
  required_cols <- c("chr", "pos", "LOD")
  if (!all(required_cols %in% names(results_df))) {
    warning("Annotation skipped: Results data frame missing required columns (chr, pos, LOD).", call. = FALSE)
    return(NULL)
  }
  if (is.null(common_chroms) || length(common_chroms) == 0) {
    warning("Annotation skipped: Missing or empty common_chroms vector.", call. = FALSE)
    return(NULL)
  }

  # Add 'end' column
  first_chrom_pos <- results_df$pos[results_df$chr == results_df$chr[1]]
  if (length(first_chrom_pos) > 1) {
    inferred_res <- stats::median(diff(first_chrom_pos), na.rm = TRUE)
  } else {
    inferred_res <- NA
  }
  if(is.na(inferred_res) || inferred_res <= 0) inferred_res <- 100
  results_df$end <- results_df$pos + inferred_res - 1

  # Bins are already filtered by the calling function based on q-value or nperm=0
  bins_to_annotate_df <- results_df

  if (nrow(bins_to_annotate_df) == 0) {
    # This message should ideally be handled by the caller, but double-check
    msg("No bins provided for annotation.")
    return(NULL)
  }

  # Create GRanges object for bins to annotate
  seqlevels_style_target <- tryCatch(GenomeInfoDb::seqlevelsStyle(GenomeInfoDb::genes(txdb, columns="tx_id"))[1], error=function(e) "UCSC")
  seqlevels_style_current <- tryCatch(GenomeInfoDb::seqlevelsStyle(bins_to_annotate_df$chr)[1], error=function(e) "UCSC")

  if(seqlevels_style_current != seqlevels_style_target) {
    msg(paste("Attempting to match sequence level styles. Current:", seqlevels_style_current, "Target:", seqlevels_style_target))
    new_seqnames <- tryCatch({
      current_map <- stats::setNames(GenomeInfoDb::seqlevels(bins_to_annotate_df$chr), GenomeInfoDb::seqlevels(bins_to_annotate_df$chr))
      GenomeInfoDb::mapSeqlevels(current_map, style = seqlevels_style_target)[bins_to_annotate_df$chr]
    }, error = function(e) {
      warning("Failed to map sequence names to TxDb style: ", e$message, ". Annotation might fail.", call. = FALSE)
      return(NULL)
    })
    if (!is.null(new_seqnames) && !anyNA(new_seqnames)) {
      bins_to_annotate_df$chr <- new_seqnames
      msg("Successfully mapped sequence names.")
    } else {
      warning("Proceeding with original sequence names despite potential style mismatch.", call.=FALSE)
    }
  }

  gr_bins <- tryCatch({
    GenomicRanges::makeGRangesFromDataFrame(
      bins_to_annotate_df,
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqnames.field = "chr",
      start.field = "pos",
      end.field = "end"
    )
  }, error = function(e) {
    warning("Failed to create GRanges object for bins: ", e$message, call. = FALSE)
    return(NULL)
  })

  if (is.null(gr_bins)) return(NULL)

  # Get gene locations from TxDb
  genes_gr <- tryCatch({
    GenomicFeatures::genes(txdb)
  }, error = function(e) {
    warning("Failed to retrieve genes from TxDb object: ", e$message, call. = FALSE)
    return(NULL)
  })

  if (is.null(genes_gr)) return(NULL)

  # Ensure seqlevels match
  common_gr_seqlevels <- intersect(GenomeInfoDb::seqlevels(gr_bins), GenomeInfoDb::seqlevels(genes_gr))
  if (length(common_gr_seqlevels) == 0) {
    warning("No common sequence levels (chromosomes) between bins and TxDb genes. Check chromosome naming.", call. = FALSE)
    return(NULL)
  }
  suppressWarnings({
    GenomeInfoDb::seqlevels(gr_bins, pruning.mode="coarse") <- common_gr_seqlevels
    GenomeInfoDb::seqlevels(genes_gr, pruning.mode="coarse") <- common_gr_seqlevels
  })

  # Find overlaps
  ov <- tryCatch({
    GenomicRanges::findOverlaps(gr_bins, genes_gr, ignore.strand = TRUE)
  }, error = function(e) {
    warning("Error finding overlaps between bins and genes: ", e$message, call. = FALSE)
    return(NULL)
  })

  if (is.null(ov) || length(ov) == 0) {
    msg("No overlapping genes found for selected bins.")
    return(NULL)
  }

  # Get gene IDs
  if ("gene_id" %in% colnames(S4Vectors::mcols(genes_gr))) {
    overlapping_gene_ids <- S4Vectors::mcols(genes_gr)$gene_id[S4Vectors::subjectHits(ov)]
  } else {
    overlapping_gene_ids <- names(genes_gr)[S4Vectors::subjectHits(ov)]
  }
  unique_gene_ids <- unique(overlapping_gene_ids)
  unique_gene_ids <- unique_gene_ids[!is.na(unique_gene_ids)]

  if (length(unique_gene_ids) == 0) {
    msg("No valid unique gene IDs found after overlap.")
    return(NULL)
  }

  # Get gene information
  key_type <- "ORF"
  if (all(grepl("^Y[A-P][LR]\\d{3}[CW](-[A-Z])?$", unique_gene_ids))) {
    key_type <- "ORF"
    msg("Using 'ORF' as keytype for annotation.")
  } else if (all(grepl("^ENS", unique_gene_ids))) {
    key_type <- "ENSEMBL"
    msg("Using 'ENSEMBL' as keytype for annotation.")
  } else {
    msg("Gene ID format unclear, attempting annotation with default keytype '", key_type, "'. If this fails, check gene ID format in your TxDb.")
  }

  gene_info <- tryCatch({
    AnnotationDbi::select(
      orgdb,
      keys = unique_gene_ids,
      columns = c("GENENAME", "DESCRIPTION"),
      keytype = key_type
    )
  }, error = function(e) {
    warning("Annotation failed using keytype '", key_type, "': ", e$message,
            "\nAvailablekeytypes might include: ", paste(AnnotationDbi::keytypes(orgdb), collapse=", "), call. = FALSE)
    if (key_type == "ORF" && !all(grepl("^Y[A-P][LR]", unique_gene_ids))) {
      msg("Retrying annotation with keytype 'GENENAME'.")
      key_type <- "GENENAME"
      tryCatch({
        AnnotationDbi::select(
          orgdb, keys = unique_gene_ids,
          columns = c("GENENAME", "DESCRIPTION"), keytype = key_type
        )
      }, error = function(e2) {
        warning("Annotation also failed using keytype 'GENENAME': ", e2$message, call. = FALSE)
        return(NULL)
      })
    } else {
      return(NULL)
    }
  })

  if (is.null(gene_info) || nrow(gene_info) == 0) {
    msg("Could not retrieve gene information from OrgDb for overlapping IDs.")
    return(NULL)
  }

  # Add position information back
  overlapping_bins <- gr_bins[S4Vectors::queryHits(ov)]
  if ("gene_id" %in% colnames(S4Vectors::mcols(genes_gr))) {
    gene_ids_overlapping <- S4Vectors::mcols(genes_gr)$gene_id[S4Vectors::subjectHits(ov)]
  } else {
    gene_ids_overlapping <- names(genes_gr)[S4Vectors::subjectHits(ov)]
  }

  # pull q-values only if the column is present (nperm > 0 case)
  q_col <- if ("q_value" %in% names(bins_to_annotate_df)) {
    bins_to_annotate_df$q_value[S4Vectors::queryHits(ov)]
  } else {
    rep(NA_real_, length(S4Vectors::queryHits(ov)))
  }

  overlap_map <- data.frame(
    gene_id      = gene_ids_overlapping,
    chr          = as.character(GenomeInfoDb::seqnames(overlapping_bins)),
    peak_bin_pos = BiocGenerics::start(overlapping_bins),
    peak_bin_lod = overlapping_bins$LOD,
    q_value      = q_col,
    stringsAsFactors = FALSE
  )

  overlap_map <- overlap_map[!is.na(overlap_map$gene_id), ]

  if (nrow(overlap_map) == 0) {
    msg("No valid gene overlaps remained after removing NA IDs.")
    return(NULL)
  }

  # ---- Collapse overlaps: pick max-LOD bin AND track min q-value per gene ----
  split_ov <- split(overlap_map, overlap_map$gene_id)

  best_bin_per_gene <- lapply(split_ov, function(df_gene) {
    max_lod_row <- df_gene[which.max(df_gene$peak_bin_lod), , drop = FALSE]
    min_q       <- min(df_gene$q_value, na.rm = TRUE)
    max_lod_row$min_q_value <- if (is.finite(min_q)) min_q else NA_real_
    max_lod_row
  })
  best_bin_per_gene <- do.call(rbind, best_bin_per_gene)
  rownames(best_bin_per_gene) <- NULL

  # ---- Merge with OrgDb annotation ----
  gene_info_merged <- merge(
    gene_info,
    best_bin_per_gene,
    by.x = key_type, by.y = "gene_id",
    all.x = TRUE
  )
  gene_info_merged <- gene_info_merged[!is.na(gene_info_merged$chr), ]

  if (nrow(gene_info_merged) == 0) {
    msg("Merging gene info with peak bin info failed.")
    return(NULL)
  }


  # Order and return unique results
  gene_info_merged <- gene_info_merged[order(match(gene_info_merged$chr, common_chroms), gene_info_merged$peak_bin_pos), ]
  gene_info_final <- unique(gene_info_merged)

  return(gene_info_final)
}


# ---- Internal Genome-wide Plotting ----

#' Plot genome-wide multipool results using ggplot2 (Internal)
#'
#' Creates a multi-panel plot showing allele frequencies and LOD scores across
#' chromosomes, with optional gene annotations and an FDR threshold line.
#' If permutations were run, genes meeting the q-value threshold are candidates
#' for labeling. If no permutations, genes near the highest LOD peaks are candidates.
#' Labels the top N candidate genes per chromosome.
#'
#' @param results_df Data frame from `multipool` containing `chr`, `pos`,
#'   `freq1_obs`, `freq2_obs`, `freq1_fit`, `freq2_fit`, `LOD`, and potentially `q_value`.
#' @param thresholds Deprecated argument, not used.
#' @param fdr_lod_threshold Numeric LOD score threshold corresponding to the desired
#'                            q-value cutoff (or NULL/Inf if none met).
#' @param genes_annotated Data frame of annotated genes (output of `internal_annotate_qtl`).
#'                       If nperm=0, this contains all overlapping genes.
#' @param nperm Number of permutations performed. Used to decide labeling strategy.
#' @param q_threshold The q-value threshold used for significance if nperm > 0.
#' @param mode Analysis mode ("replicates" or "contrast") for axis labeling.
#' @param chrom_order Character vector specifying the desired order of chromosomes for plotting.
#' @param n_label Integer, maximum number of genes to label per chromosome.
#' @return A ggplot object or NULL if plotting fails.
#' @keywords internal
internal_plot_genome_wide <- function(results_df, thresholds = NULL, fdr_lod_threshold = NULL, genes_annotated, nperm, q_threshold, mode, chrom_order, n_label = N_GENES_LABEL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed, cannot plot. Install with install.packages('ggplot2')", call. = FALSE)
    return(NULL)
  }
  use_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
  if (!use_ggrepel) {
    warning("ggrepel not installed, gene labels may overlap. Install with install.packages('ggrepel')", call. = FALSE)
  }

  if (nrow(results_df) == 0) {
    warning("Plotting skipped: Input results data frame is empty.", call. = FALSE)
    return(NULL)
  }
  required_cols <- c("chr", "pos", "freq1_obs", "freq2_obs", "freq1_fit", "freq2_fit", "LOD")
  if (!all(required_cols %in% names(results_df))) {
    warning("Plotting skipped: Results data frame missing required columns.", call. = FALSE)
    return(NULL)
  }

  # --- Data Preparation ---
  results_df$pos <- as.numeric(results_df$pos)
  results_df$freq1_obs <- as.numeric(results_df$freq1_obs)
  results_df$freq2_obs <- as.numeric(results_df$freq2_obs)
  results_df$freq1_fit <- as.numeric(results_df$freq1_fit)
  results_df$freq2_fit <- as.numeric(results_df$freq2_fit)
  results_df$LOD <- as.numeric(results_df$LOD)
  plot_df <- results_df[!is.na(results_df$pos) & !is.na(results_df$LOD), ]
  if (nrow(plot_df) == 0) {
    warning("Plotting skipped: No valid data points after removing NAs.", call. = FALSE)
    return(NULL)
  }
  plot_chroms <- unique(plot_df$chr)
  missing_chroms <- setdiff(plot_chroms, chrom_order)
  if (length(missing_chroms) > 0) {
    warning("Chromosomes found in results but not in provided chrom_order: ",
            paste(missing_chroms, collapse=", "), ". Appending them to the order.", call. = FALSE)
    chrom_order <- c(chrom_order, missing_chroms)
  }
  chrom_order <- intersect(chrom_order, plot_chroms)
  if (length(chrom_order) == 0) {
    warning("Plotting skipped: No common chromosomes between provided order and results.", call. = FALSE)
    return(NULL)
  }
  plot_df$chr <- factor(plot_df$chr, levels = chrom_order)

  # --- Scaling ---
  lod_values_for_max <- c(0, plot_df$LOD, if (!is.null(fdr_lod_threshold) && is.finite(fdr_lod_threshold)) fdr_lod_threshold else numeric(0))
  max_lod <- max(lod_values_for_max[is.finite(lod_values_for_max)], na.rm = TRUE)
  if (!is.finite(max_lod) || max_lod < 1e-9) max_lod <- 1
  # If threshold is very high, ensure plot goes slightly above it
  if (!is.null(fdr_lod_threshold) && is.finite(fdr_lod_threshold) && fdr_lod_threshold > max_lod) {
    max_lod <- fdr_lod_threshold * 1.1
  }
  freq_values_for_max <- c(0, 1, plot_df$freq1_obs, plot_df$freq2_obs,
                           plot_df$freq1_fit, plot_df$freq2_fit)
  max_freq <- max(freq_values_for_max[is.finite(freq_values_for_max)], na.rm = TRUE)
  max_freq <- max(1.0, max_freq * 1.05)
  scale_factor <- max_freq / max_lod
  plot_df$LOD_scaled <- plot_df$LOD * scale_factor

  # --- Base Plot ---
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pos / 1e6)) +
    ggplot2::geom_point(ggplot2::aes(y = freq1_obs), color = "lightblue", alpha = 0.4, shape = 1, size = 1.5,
                        na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(y = freq2_obs), color = "lightcoral", alpha = 0.4, shape = 1, size = 1.5,
                        na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = freq1_fit), color = "blue", linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = freq2_fit), color = "red", linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = LOD_scaled), color = "black", linewidth = 1.0) +
    ggplot2::facet_wrap(~ chr, scales = "free_x", ncol = ceiling(length(chrom_order)/4)) +
    ggplot2::scale_y_continuous(
      name = "Allele Frequency", limits = c(0, max_freq),
      expand = ggplot2::expansion(mult = c(0, 0.05)),
      sec.axis = ggplot2::sec_axis(~ . / scale_factor, name = "LOD Score")
    ) +
    ggplot2::scale_x_continuous(name = "Position (Mb)", expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::ggtitle(paste("Multipool Genome Scan Results - Mode:", mode)) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(linewidth = 0.3, colour = "grey90"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey85", color = "grey50"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size=8),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )

  # --- Add FDR Threshold Line ---
  if (!is.null(fdr_lod_threshold) && is.finite(fdr_lod_threshold)) {
    threshold_scaled <- fdr_lod_threshold * scale_factor
    # Only plot if it's within the y-axis range
    if(threshold_scaled <= max_freq && threshold_scaled >= 0) {
      p <- p + ggplot2::geom_hline(
        yintercept = threshold_scaled,
        color = "darkorange", # Use a different color for FDR
        linetype = "dashed",
        linewidth = 0.8
      ) +
        ggplot2::annotate( # Use annotate for a single label across facets
          geom = "text", x = Inf, y = threshold_scaled,
          label = paste0("q < ", q_threshold),
          hjust = 1.05, vjust = -0.5, color = "darkorange", size = 2.5
        )
    }
  }

  # --- Add Gene Annotations --------------------------------------------------
  if (!is.null(genes_annotated) && nrow(genes_annotated) > 0) {

    req_cols <- c("chr", "peak_bin_pos", "GENENAME", "peak_bin_lod")
    if (!all(req_cols %in% names(genes_annotated))) {
      warning("Gene annotation plotting skipped – missing ",
              paste(setdiff(req_cols, names(genes_annotated)), collapse = ", "),
              call. = FALSE)
    } else {

      ## ------------------------------------------------------------------
      ## guarantee a gene-level q_value column exists BEFORE we go further
      ## ------------------------------------------------------------------
      if ("q_value" %in% names(genes_annotated)) {
        # fine
      } else if ("min_q_value" %in% names(genes_annotated)) {
        genes_annotated$q_value <- genes_annotated$min_q_value
      } else {
        genes_annotated$q_value <- NA_real_
      }

      genes_to_label <- NULL

      ## 1)  PERMUTATIONS  -------------------------------------------------
      if (nperm > 0) {

        res_cols <- c("chr", "pos")
        if ("q_value" %in% names(results_df))
          res_cols <- c(res_cols, "q_value")

        gmerged <- merge(
          genes_annotated,
          results_df[, res_cols, drop = FALSE],
          by.x = c("chr", "peak_bin_pos"),
          by.y = c("chr", "pos"),
          all.x = TRUE
        )

        ## ------ normalise any suffixed names from merge() ----------------
        if (!"q_value" %in% names(gmerged)) {
          if ("q_value.x" %in% names(gmerged))
            gmerged$q_value <- gmerged$q_value.x
          else if ("q_value.y" %in% names(gmerged))
            gmerged$q_value <- gmerged$q_value.y
          else
            gmerged$q_value <- NA_real_
        }

        sig <- subset(gmerged, !is.na(q_value) & q_value <= q_threshold)

        if (nrow(sig) > 0) {
          genes_to_label <- do.call(
            rbind,
            lapply(split(sig, sig$chr), function(d) {
              d <- d[order(d$peak_bin_lod, decreasing = TRUE, na.last = NA), ]
              utils::head(d, n_label)
            })
          )
          rownames(genes_to_label) <- NULL
        }

        ## 2)  NO PERMUTATIONS  ---------------------------------------------
      } else {
        msg("nperm = 0: labeling genes near top LOD peaks.")
        genes_to_label <- do.call(
          rbind,
          lapply(split(genes_annotated, genes_annotated$chr), function(d) {
            d <- d[order(d$peak_bin_lod, decreasing = TRUE, na.last = NA), ]
            utils::head(d, n_label)
          })
        )
        if (!is.null(genes_to_label)) rownames(genes_to_label) <- NULL
      }

      ## 3)  DRAW ----------------------------------------------------------
      if (!is.null(genes_to_label) && nrow(genes_to_label) > 0) {

        gene_plot_df        <- genes_to_label
        gene_plot_df$chr    <- factor(gene_plot_df$chr,
                                      levels = levels(plot_df$chr))
        gene_plot_df$pos_mb <- gene_plot_df$peak_bin_pos / 1e6
        gene_plot_df        <- gene_plot_df[!is.na(gene_plot_df$chr), ]

        if (nrow(gene_plot_df) > 0) {
          geom_fun   <- if (use_ggrepel) ggrepel::geom_text_repel
          else             ggplot2::geom_text
          base_args  <- list(
            mapping = ggplot2::aes(x = pos_mb,
                                   y = max_freq * 0.95,
                                   label = GENENAME),
            data    = gene_plot_df,
            color   = "purple",
            size    = 2.5,
            max.overlaps = Inf
          )
          extra_args <- if (use_ggrepel)
            list(segment.color = "grey50",
                 segment.size  = 0.3,
                 min.segment.length = 0)
          else
            list(hjust = 0.5, vjust = 0)

          p <- p + do.call(geom_fun, c(base_args, extra_args, list(na.rm = TRUE)))
        }
      } else {
        msg("No genes selected for labeling (nothing to plot).")
      }
    }
  }




  return(p)
}

# ---- END OF SCRIPT ----
