% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multipool.R
\name{internal_annotate_qtl}
\alias{internal_annotate_qtl}
\title{Annotate QTL regions with yeast genes (Internal Implementation)}
\usage{
internal_annotate_qtl(
  results_df,
  q_threshold = NULL,
  txdb,
  orgdb,
  common_chroms
)
}
\arguments{
\item{results_df}{Data frame of bins to annotate (pre-filtered based on significance or all bins if nperm=0).
Must include chr, pos, LOD.}

\item{q_threshold}{Optional numeric q-value significance threshold (only used for messaging).}

\item{txdb}{TxDb object (e.g., TxDb.Scerevisiae.UCSC.sacCer3.sgdGene).}

\item{orgdb}{OrgDb object (e.g., org.Sc.sgd.db).}

\item{common_chroms}{Character vector of common chromosome names used in analysis (for ordering).}
}
\value{
Data frame of unique annotated genes, or NULL.
}
\description{
Finds genes overlapping bins. If a q_threshold is provided (implicitly via
results_df containing a q_value column), only genes overlapping bins
meeting the threshold are returned. Otherwise (if nperm=0), all
overlapping genes are returned.
}
\keyword{internal}
