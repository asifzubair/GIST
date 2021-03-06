% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{preprocess_expr_mat}
\alias{preprocess_expr_mat}
\title{Impute, filter and normalize expression matrix}
\usage{
preprocess_expr_mat(
  counts_mat,
  impute = TRUE,
  k = 5,
  d = 10,
  seed = 42,
  doQC = FALSE,
  normalize = NULL,
  only_hvg = TRUE
)
}
\arguments{
\item{counts_mat}{Dataframe of counts, cells are columns, genes are rows}

\item{impute}{Logical indicating whether imputation should be done (default: TRUE)}

\item{k}{Number of nearest neighbours to aggregate (default: 5)}

\item{d}{Number of principal components (default: 10)}

\item{seed}{seed to be passed to knnSmooth algorithm}

\item{doQC}{Should we do filtering (default: FALSE)}

\item{normalize}{Do one of 'scale, scale-quantile, quantile, sct' normalization (default: None)}

\item{only_hvg}{Return only HVG from SCTransform (default: TRUE)}
}
\value{
a matrix of gene counts
}
\description{
Given a matrix of SC/ST counts, this function will
impute gene counts using knn_smoothing
filter based on criteria that gene has
2 reads in atleast 10 cells
and finally, normalize the matrix using the chosen method
}
