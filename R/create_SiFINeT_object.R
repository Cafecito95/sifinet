#' The SiFINeT Class
#'
#' @slot data a list of cell (row) by gene (column) count matrices ("counts" and "counts_subcohort")
#' @slot sparse whether the count matrices should be analyzed as sparse matrices
#' @slot meta.data matrix of meta data, the number of rows should equal to the number of cells
#' @slot gene.name a vector of names of genes with length equal to the number of genes 
#' @slot data.name name of the dataset
#' @slot data.name_subcohort name of the subcohort dataset
#' @slot n number of cells in the dataset
#' @slot n_subcohort number of cells in the subcohort dataset
#' @slot p number of genes in the dataset
#' @slot p_subcohort number of genes in the subcohort dataset
#' @slot data.thres binarized count matrix 
#' @slot coexp matrix of genes coexpression
#' @slot est_ms estimated mean and sd of coexpression values
#' @slot thres lower bound of coexpression (or absolute value of coexpression) for network edge assignment
#' @slot q5 50% quantile for each gene
#' @slot kset index of kept genes after the filtering step
#' @slot conn list of connectivities in absolute network
#' @slot conn2 list of connectivities in positive sub-network
#' @slot fg_id index of the candidate feature genes
#' @slot uni_fg_id index of the candidate unique feature genes
#' @slot uni_cluster cluster result of the candidate unique feature genes
#' @slot selected_cluster selected unique feature gene clusters
#' @slot featureset detected set of feature genes
#' 
#' @name SiFINeT-class
#' @exportClass SiFINeT
#'
SiFINeT <- setClass(
  Class = 'SiFINeT',
  slots = c(
    data = 'list',
    sparse = 'logical',
    meta.data = 'matrix',
    gene.name = 'vector',
    data.name = 'character',
    data.name_subcohort = 'character',
    n = 'numeric',
    n_subcohort = 'numeric',
    p = 'numeric',
    p_subcohort = 'numeric',
    data.thres = 'list',
    coexp = 'matrix',
    est_ms = 'list',
    thres = 'numeric',
    q5 = 'numeric',
    kset = 'integer',
    conn = 'list',
    conn2 = 'list',
    fg_id = 'integer',
    uni_fg_id = 'integer',
    uni_cluster = 'numeric',
    selected_cluster = 'numeric',
    featureset = 'list'
  )
)

#' create_SiFINeT_object
#' 
#' The function classifies count data based on thresholds 
#' defined by quantile regression
#' @param counts primary count matrix
#' @param counts_subcohort sub-cohort count matrix
#' @param gene.name name of the features
#' @param meta.data data.frame of meta data
#' @param data.name name of dataset
#' @param sparse whether the count matrices should be analyzed as sparse matrices
#' @param rowfeature whether the count matrices are feature (row) by cell (column) 
#' @return a SiFINeT object
#' @export
#' 
create_SiFINeT_object <- function(counts, counts_subcohort, gene.name = NULL, 
                                  meta.data = NULL, data.name = NULL, data.name_subcohort = NULL,
                                  sparse = FALSE, rowfeature = TRUE){
  if (rowfeature){
    counts <- t(counts)
    counts_subcohort <- t(counts_subcohort)
  }
  if (is.null(gene.name)){
    gene.name <- colnames(counts)
  }
  if (is.null(data.name)){
    data.name <- "data1"
    data.name_subcohort <- "data2"
  }
  if (is.null(meta.data)){
    meta.data <- matrix(0, nrow(counts), 0)
  }
  data <- list(counts = counts, counts_subcohort = counts_subcohort)
  names(data) <- c(data.name, data.name_subcohort)
  object <- new(
    Class = 'SiFINeT',
    data = data,
    sparse = sparse,
    meta.data = meta.data,
    gene.name = gene.name,
    data.name = data.name,
    data.name_subcohort = data.name_subcohort,
    n = nrow(counts),
    n_subcohort = nrow(counts_subcohort),
    p = ncol(counts),
    p_subcohort = ncol(counts_subcohort),
    q5 = apply(counts, 2, quantile, 0.5),
    kset = 1:ncol(counts),
    featureset = list(unique = list(),
                      shared = list(),
                      enriched = list())
  )
  return(object)
}
