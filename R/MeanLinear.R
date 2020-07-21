#' @title
#' Linear Mean Function
#'
#' @description
#' A linear mean.
#'
#' @details
#' This function gives a linear mean; denoting its hyper vector by \eqn{\beta},
#' for an input matrix \eqn{X}, the prior mean is \eqn{X \beta}.
#'
#' @export
MeanLinear <- R6::R6Class(
    classname = "MeanLinear",
    inherit = MeanFunction,
    public = list(

        ## Data members
        #' @field name A character vector of length one giving the mean
        #'     function's name; "linear"
        name = "linear",
        #' @field hypers A numeric vector giving the mean function's
        #'     hyperparameters
        hypers = numeric(),

        ## Methods
        #' @description
        #' Compute function prior mean
        #' @param X The input values (should be a numeric matrix)
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     mean function. If NULL (the default), the hypers data
        #'     member is used.
        mean = function(X, hypers = NULL) {
            if ( is.null(hypers) ) {
                hypers <- self$hypers
            }
            return(X %*% hypers)
        },
        #' @description
        #' Compute partial derivatives of mean function with respect to
        #' its hyperparameters
        #' @param X The nput values (should be a numeric matrix)
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     mean function. If NULL (the default), the hypers data
        #'     member is used.
        #' @param param An integer vector of length one; which element of
        #'     \code{hypers} should the derivative be taken with respect to?
        #'     The default is 1
        parameter_derivative = function(X, hypers = NULL, param = 1) {
            if ( param > ncol(X) ) {
                stop("param must be <= ncol(X).")
            }
            return(X[ , param])
        },
        #' @description
        #' Compute partial derivatives of mean function with respect to
        #' its inputs
        #' @param X The input values (should be a numeric matrix)
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     mean function. If NULL (the default), the hypers data
        #'     member is used.
        #' @param dimension an integer vector of length one giving the dimension
        #'     of X with respect to which the derivative is being taken; the
        #'     default is 1
        input_derivative = function(X, hypers = NULL, dimension = 1) {
            if ( dimension > ncol(X) ) {
                stop("dimension must be <= ncol(X).")
            }
            if ( is.null(hypers) ) {
                hypers <- self$hypers
            }
            return(rep(hypers[dimension], nrow(X)))
        },

        ## Constructor
        #' @description
        #' Create a new MeanLinear object
        #' @param hypers A numeric vector giving hyperparameters for
        #'     the mean function.
        initialize = function(hypers = numeric()) {
            self$hypers = hypers
        }
    )
)
