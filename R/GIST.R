#' Guiding-Image Spatial Transcriptomics
#'
#' @param st_expression Spatial transcriptomics expression matrix (genes by spot)
#' @param sig_mat Signature matrix (genes by cell type)
#' @param prior_values Vector of prior values
#' @param prior_index Index in signature matrix of cell type which needs a prior on it
#' @param prior_lambda Tunable hyperparameter \eqn{\lambda} for the prior beta distribution
#' @param num_cores How many cores to use for parallel processing
#' @param ... Arguments to be passed to \code{rstan::sampling} (e.g. \code{chains, iter, init, verbose, refresh})
#' @return An object of type dataframe with posterior mean estimates
#' @export
GIST <- function(st_expression, sig_mat, prior_values = NULL, prior_index = NULL,  prior_lambda = 50, num_cores = 1, ...){

  # TODO: this following sanity check is very basic,
  # TODO: need to handle situations in which either row or column vectors are passed
  if (is.numeric(st_expression)) st_expression <- as.data.frame(st_expression)

  genes <- intersect(rownames(sig_mat), rownames(st_expression))
  message('Overalp betweeb ST and SC genes:', length(genes), '\n')

  st_expression <- st_expression[genes, ]
  sig_mat <- sig_mat[genes, ]
  stopifnot(all(rownames(sig_mat) == rownames(st_expression)))

  numGenes = nrow(st_expression)
  numCellTypes = ncol(sig_mat)

  #TODO: ideally, I want to do some diagnostic checks based on Stan warnings.
  #TODO: will probably write a curried function that does diagnostic checks.
  if (is.null(prior_index))
    model = stanmodels$GIST_base_model
  else
    model = stanmodels$GIST_enhanced_model

  ## get kwargs and register parallel backend
  args0 = list(...)
  if (num_cores > 1){
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
  }

  if(foreach::getDoParRegistered())
    message("Parallel execution would follow metrics:",
            "\nName --> ", foreach::getDoParName(),
            "\nVersion --> ", foreach::getDoParVersion(),
            "\nWorkers --> ", foreach::getDoParWorkers(), "\n")
  else
    message("No parallel backend detected")

  parout <- foreach::foreach(i = seq_len(NCOL(st_expression)),
                             .packages = "rstan",
                             .combine = rbind) %dopar%  {
      if (is.null(prior_index))
        standata <- list(numGenes = numGenes, numCellTypes = numCellTypes,
                         exprMixVec = st_expression[, i], sigMat = sig_mat)
      else
        standata <- list(numGenes = numGenes, numCellTypes = numCellTypes,
                         exprMixVec = st_expression[, i], sigMat = sig_mat,
                         prior_phi = prior_values[i], prior_index = prior_index, prior_lambda = prior_lambda)
      print(model)
      nmfOut <- do.call(sampling, c(list(object = model, data = standata), args0))
      stanSumNmf <- rstan::extract(nmfOut)$estimatedProportionsVecSimp[, 1:numCellTypes]
      pEstimatesList <- colMeans(stanSumNmf)

      return(pEstimatesList)
    }

  if (num_cores > 1) parallel::stopCluster(cl)

  colnames(parout) <- colnames(sig_mat)
  rownames(parout) <- colnames(st_expression)
  return(parout)
}
