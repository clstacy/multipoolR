% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{multipool_example_data}
\alias{multipool_example_data}
\alias{multipoolR_example_pool1}
\alias{multipoolR_example_pool2}
\title{Example Multipool Data with Sloped QTLs}
\format{
A list containing two data frames, which are loaded into the
environment when \code{data(multipool_example_data)} is called:
\describe{
\item{multipoolR_example_pool1}{Data frame for Pool 1 (7000 rows)}
\item{multipoolR_example_pool2}{Data frame for Pool 2 (7000 rows)}
}
Each data frame has columns:
\describe{
\item{chr}{Chromosome name (chrIII or chrV)}
\item{pos}{Genomic position (integer)}
\item{a}{Allele count for reference/parent 1 (integer)}
\item{b}{Allele count for reference/parent 2 (integer)}
}
}
\source{
Simulated data generated for package examples. See README for generation code.
}
\usage{
data(multipool_example_data)
}
\description{
A simulated dataset containing allele counts for two pools (1 and 2)
across two chromosomes (chrIII and chrV), designed to illustrate
the usage of the multipoolR package.
}
\details{
This dataset includes simulated QTLs with effects that ramp linearly
from the baseline frequency (0.5) at the edges to a peak frequency
at the center of the defined QTL region.

\strong{Pool 1 QTLs:}
\itemize{
\item chrIII: Peak frequency 0.80 between 140kb and 175kb.
\item chrV: Peak frequency 0.35 between 400kb and 420kb.
}

\strong{Pool 2 QTLs:}
\itemize{
\item chrIII: Peak frequency 0.65 between 140kb and 175kb.
\item chrV: Peak frequency 0.60 between 400kb and 420kb.
}

The data was generated using the \code{generate_multipool_data} function (see README
or examples for its definition) with specific seeds for reproducibility.
}
\keyword{datasets}
