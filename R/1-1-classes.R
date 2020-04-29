#################################
#### Psiform class
#' An S4 class containing predictors (x), response (y) and sample weights (w)
#'
#' @slot cluster A numeric vector of cluster number
#' @slot param A list of parameters in the RZiMM model
#' @slot importance A numeric vector for the importance of each gene
#' @slot info A list of items, including log-likelihood, bic, and etc
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
