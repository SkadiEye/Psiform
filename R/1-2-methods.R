###########################################################
### Define generic functions: show, print, '[', '$', plot

#################################
#### Psiform class

#' @name show
#' @rdname show
#'
#' @title Method show for the package
#'
#' @description The method show for \code{Psiform} object.
#'
#' @param object A \code{Psiform} object.
#'
#' @seealso
#' \code{\link{Psiform-class}}\cr
#'
NULL

#' @rdname show
#' @importFrom methods show
#'
#' @export
setMethod("show",
          "Psiform",
          function(object) {

            cat("c:", object@c, "Object \n")
            cat("beta:", object@beta[1:min(6, length(object@beta))], "(first 6 cells) \n")
          })

