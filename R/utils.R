#' read in MEX format data produced by spaceranger
#' source:
#' https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/matrices
#' @param matrix_dir directory containing counts, barcode and features
load10xMatrix <- function(matrix_dir) {
  barcode.path <- file.path(matrix_dir, 'barcodes.tsv')
  features.path <- file.path(matrix_dir, 'features.tsv')
  matrix.path <- file.path(matrix_dir, 'matrix.mtx')
  mat <- Matrix::readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V2
  return(mat)
}

#' read in the tissue positions file produced by spaceranger
#' ensures consistent columnn nomenclature
#' @param tp_file tissue positions fils from spaceranger
load10xTP <- function(tp_file) {
  columns <- c('barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_col', 'pxl_row')
  tp <- read.csv(tp_file,
                 col.names = columns,
                 row.names = 1)
  return(tp)
}

#' Impute, filter and normalize
#' Given a matrix of SC/ST counts, this will
#' impute gene counts using knn_smoothing
#' filter based on criteria that gene has
#' 2 reads in atleast 10 cells
#' and finally, quantile normalize the matrix
#' @param counts_mat a dataframe of counts, cells are columns, genes are rows
#' @param impute
#' @param k
#' @param d
#' @param doQC do some filtering ?
#' @param normalize do one of 'scale, scale-quantile, quantile, sct' normalization (default: None)
#' @param only_hvg return only HVG from SCTransform (default: TRUE)
preprocess_expr_mat <- function(counts_mat, impute = T, k = 5, d = 10, seed = 42,
                                doQC = F, normalize = NULL, only_hvg = T){

  if (impute){
    if(!exists('knn_smoothing'))
      stop('make sure you have sourced the imputation file')
    message('Smoothing with ', k, ' neighbours and ', d, ' dimensions ...')
    counts_mat <- knn_smoothing(as.matrix(counts_mat),
                                k = k, d = d,
                                seed = seed)
  }

  if (doQC) {
    message('performing filtering ...')
    genes <- which(rowSums(counts_mat > 2) > 10)
    counts_mat <- counts_mat[genes, ]
  }

  if (!is.null(normalize)){
    stopifnot(normalize %in% c('scale', 'quantile', 'scale-quantile', 'sct'))
    message('performing ', normalize, ' normalization ...')
    counts_mat <- normalize.decon(counts_mat,
                                  doQC = F,
                                  method = normalize, only_hvg = only_hvg)
  }

  return(counts_mat)

}


#' normalize a genes by cell matrix
#' optionally perfrom some QC
#' @param data a genes x cell matrix
#' @param doQC should QC be performed: min.cells = 10%, min.features = 500
#' @param method choose a method from 'scale', 'quantile', 'scale-quantile', 'sct'
#' @param only_hvg if using method 'sct', return only highly variable genes after normalization
normalize.decon <- function(data,
                            doQC = FALSE,
                            method = c('scale', 'quantile', 'scale-quantile', 'sct'),
                            only_hvg = T){

  method = match.arg(method)
  if (method %in% c('scale', 'scale-quantile'))
    scale.factor = median(colSums(data))

  if (doQC){
    data <- Seurat::CreateSeuratObject(data,
                                       min.cells = 0.1*ncol(data),
                                       min.features = 500,
                                       assay = 'RNA')
  } else {
    data <- Seurat::CreateSeuratObject(data,
                                       min.cells = 0,
                                       min.features = 0,
                                       assay = 'RNA')
  }

  if (method == 'scale'){
    data <- Seurat::NormalizeData(data,
                                  normalization.method = 'RC',
                                  scale.factor = scale.factor)
    data <- as.data.frame(Seurat::GetAssayData(data,
                                               assay = 'RNA',
                                               slot = 'data'))
  } else if (method == 'quantile'){
    data <- as.data.frame(Seurat::GetAssayData(data,
                                               assay = 'RNA',
                                               slot = 'counts'))
    data <- normalize.quantiles2(data)

  } else if (method == 'scale-quantile'){
    data <- Seurat::NormalizeData(data,
                                  normalization.method = 'RC',
                                  scale.factor = scale.factor)
    data <- as.data.frame(Seurat::GetAssayData(data,
                                               assay = 'RNA',
                                               slot = 'data'))
    data <- normalize.quantiles2(data)
  } else if (method == 'sct') {
    if (only_hvg)
      message('SCT will return only HVG ...')
    data <- Seurat::SCTransform(data, assay = 'RNA', return.only.var.genes = only_hvg)
    data <- as.data.frame(Seurat::GetAssayData(data,
                                               assay = 'SCT',
                                               slot = 'scale.data'))
  }

  return(data)

}


#' Normalize quantiles wrapper
#' to handle loss of rownames
#' and colnames
#' @param counts_mat a matrix of couts to be normalized
normalize.quantiles2 <- function(counts_mat){
  RNs <- rownames(counts_mat)
  CNs <- colnames(counts_mat)
  convert = F
  if (!is.matrix(counts_mat)) {
    convert = T
    counts_mat <- as.matrix(counts_mat)
  }
  counts_mat <- preprocessCore::normalize.quantiles(counts_mat)
  rownames(counts_mat) = RNs
  colnames(counts_mat) = CNs
  if (convert) counts_mat <- as.data.frame(counts_mat)
  return(counts_mat)
}


#' make signature matrix from
#' imputed, filtered, normalized
#' single cell counts
#' @param sc_counts a matrix/dataframe of single cell counts, with cells as columns and genes as rows
#' @param sc_labels a matrix/dataframe of cell type labels, two cols, one for cell and one for cell type
#' @param save_file optionally save the signature matrix
make_signature_matrix <- function(sc_counts, sc_labels, compute_std = F, save_prefix = NULL){

  if (is.matrix(sc_counts)) sc_counts <- as.data.frame(sc_counts)
  if (is.matrix(sc_labels)) sc_labels <- as.data.frame(sc_labels)
  colnames(sc_labels) <- c('cell', 'bio_celltype')

  sig_mat <- sc_counts %>%
    t %>% as.data.frame %>%
    tibble::rownames_to_column('cell') %>%
    dplyr::inner_join(sc_labels, by = 'cell') %>%
    dplyr::select(bio_celltype, everything()) %>%
    dplyr::select(-cell)

  sig_mat_means <- sig_mat %>%
    dplyr::group_by(bio_celltype) %>%
    dplyr::summarise_all(mean) %>%
    tibble::column_to_rownames('bio_celltype') %>%
    t %>% as.data.frame()

  sig_mat_std <- NULL
  if (compute_std){
    sig_mat_std <- sig_mat %>%
      dplyr::group_by(bio_celltype) %>%
      dplyr::summarise_all(sd) %>%
      tibble::column_to_rownames('bio_celltype') %>%
      t %>% as.data.frame()
  }

  if(!is.null(save_prefix)){
    saveRDS(sig_mat_means, file = paste0(save_prefix, '_means.rds'))
    if (compute_std)
      saveRDS(sig_mat_std, file = paste0(save_prefix, '_std.rds'))
  }

  invisible(list(sig_means = sig_mat_means, sig_std = sig_mat_std))
}


#' visualization of spots expression
#' across tissue for a single variable
#' @param df spot by value matrix, spots are rows, values are cols
#' @param tp matrix of spot locations, should contain
#'           pxl_col, pxl_row indication column and row location of pixel
#' @param trans any transofrmation that needs to be done to values - log10 etc.
#' @param legend.label label for the value legend
plot.spots0 <- function(df, tp, trans = 'identity', legend.label = 'marker exp.', size_m = 4){
  p <- df %>%
    merge_with_locations(tp) %>%
    ggplot(aes(pxl_row, -pxl_col, color = marker)) +
    geom_point(size = size_m) +
    scale_colour_gradient(low="#f8d39f", high = "#800000", trans = trans) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 14),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.background = element_rect(fill = "black")) +
    labs(colour = legend.label)
  return(p)
}


#' visualization of spots expression
#' across tissue
#' @param df spot by value matrix, spots are rows, values are cols
#' @param tp matrix of spot locations, should contain
#'           pxl_col, pxl_row indication column and row location of pixel
#' @param trans any transofrmation that needs to be done to values - log10 etc.
#' @param legend.label label for the value legend
#' @param joint_scale plot common scale for all spatial distributions
plot.spots <- function(df, tp, trans = 'identity',
                       legend.label = 'marker exp.', ncols = 3, joint_scale = FALSE, my_theme = NULL, size_m = 1){
  default_theme <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #strip.text = element_text(size = 14),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "black"))
  if (is.null(my_theme)) my_theme = default_theme
  df <- df %>%
    merge_with_locations(tp)
  if (!joint_scale){
    p <-  df %>%
      group_split(type) %>%
      purrr::map(
        ~ggplot(., aes(pxl_row, -pxl_col, color = marker)) +
          geom_point(size = size_m, shape=16) +
          scale_colour_gradient(low="#f8d39f", high = "#800000", trans = trans) +
          facet_wrap(~ type, labeller = function(x) label_value(x, multi_line = FALSE)) +
          labs(colour = legend.label) +
          my_theme
      ) %>%
      cowplot::plot_grid(plotlist = ., align = 'hv', ncol = ncols)
  } else {
    p <- df %>% ggplot(aes(pxl_row, -pxl_col, color = marker)) +
      geom_point(size = size_m, shape=16) +
      scale_colour_gradient(low="#f8d39f", high = "#800000", trans = trans) +
      facet_wrap(~ type, labeller = function(x) label_value(x, multi_line = FALSE), ncol = ncols) +
      my_theme +
      labs(colour = legend.label)
  }
  return(p)
}
