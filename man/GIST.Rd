% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GIST.R
\name{GIST}
\alias{GIST}
\title{Guiding-Image Spatial Transcriptomics}
\usage{
GIST(
  st_expression,
  sig_mat,
  prior_values = NULL,
  prior_index = NULL,
  prior_lambda = 50,
  num_cores = 1,
  ...
)
}
\arguments{
\item{st_expression}{Spatial transcriptomics expression matrix (genes by spot)}

\item{sig_mat}{Signature matrix (genes by cell type)}

\item{prior_values}{Vector of prior values}

\item{prior_index}{Index in signature matrix of cell type which needs a prior on it}

\item{prior_lambda}{Tunable hyperparameter \eqn{\lambda} for the prior beta distribution}

\item{num_cores}{How many cores to use for parallel processing}

\item{...}{Arguments to be passed to \code{rstan::sampling} (e.g. \code{chains, iter, init, verbose, refresh})}
}
\value{
An object of type dataframe with posterior mean estimates
}
\description{
Guiding-Image Spatial Transcriptomics
}
