#################################
#### Psiform class
#' An S4 class containing predictors (x), response (y) and sample weights (w)
#'
#' @slot c Phenotype values.
#' @slot beta Parameters for genes.
#' @slot u U matrix.
#' @slot v V matrix.
#' @slot sigma_sq Variances of the residual terms.
#' @slot param Model parameters.
#'
#' @export
setClass("Psiform",
         slots = list(
           c = "ANY",
           beta = "numeric",
           u = "matrix",
           v = "matrix",
           sigma_sq = "numeric",
           param = "list"
         ))
