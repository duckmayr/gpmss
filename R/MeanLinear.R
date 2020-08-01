#' @title
#' Linear Mean Function
#'
#' @description
#' A linear mean.
#'
#' @details
#' This function gives a linear mean; denoting its hyper vector by
#' \eqn{(\alpha, \beta)}, where \eqn{\alpha} is an (optional) intercept,
#' and \eqn{\beta} is a vector of coefficients whose length is the same
#' as the number of columns in \eqn{X}, then for an input matrix \eqn{X},
#' the prior mean is \eqn{\alpha + X \beta}.
#'
#' Note: Usually in regression models we would instead use the notation
#' \eqn{X \beta}, where the first element of \eqn{\beta} is a coefficient
#' for the intercept and \eqn{X} has a column of ones prepended to it.
#' The reason we do not do this in GP regression and classification is
#' because adding a one column to the design matrix is generally not
#' desirable in terms of the inputs to the kernel, and is not generally
#' helpful for other mean functions as well. So, we simply address the
#' intercept issue by adding a hyperparameter to a linear mean function
#' that is not involved with any multiplication with predictors.
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
        #' @field intercept A logical vector of length one; does the linear
        #'     mean include an intercept?
        intercept = logical(),

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
                if ( length(hypers) == 1 ) {
                    hypers <- rep(hypers, length.out = ncol(X) + self$intercept)
                    self$hypers <- hypers
                }
            }
            if ( length(hypers) != ncol(X) + self$intercept ) {
                stop("The number of coefficients provided to MeanLinear ",
                     "does not match the number of columns in X.")
            }
            if ( self$intercept ) {
                alpha <- hypers[ 1]
                beta  <- hypers[-1]
            } else {
                alpha <- 0
                beta  <- hypers
            }
            return(c(alpha + X %*% beta))
        },
        #' @description
        #' Compute partial derivatives of mean function with respect to
        #' its hyperparameters
        #' @param X The input values (should be a numeric matrix)
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     mean function. If NULL (the default), the hypers data
        #'     member is used.
        #' @param param An integer vector of length one; which element of
        #'     \code{hypers} should the derivative be taken with respect to?
        #'     The default is 1
        parameter_derivative = function(X, hypers = NULL, param = 1) {
            if ( param > (ncol(X) + self$intercept) ) {
                stop("param must be <= ncol(X).")
            }
            if ( self$intercept ) {
                if ( param == 1 ) {
                    return(rep(1, nrow(X)))
                } else {
                    param <- param - 1
                }
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
                if ( length(hypers) == 1 ) {
                    hypers <- rep(hypers, length.out = ncol(X) + self$intercept)
                    self$hypers <- hypers
                }
            }
            if ( length(hypers) != ncol(X) + self$intercept ) {
                stop("The number of coefficients provided to MeanLinear ",
                     "does not match the number of columns in X.")
            }
            return(rep(hypers[dimension + self$intercept], nrow(X)))
        },

        ## Constructor
        #' @description
        #' Create a new MeanLinear object
        #' @param hypers A numeric vector giving hyperparameters for
        #'     the mean function. If the provided hypers are of length one,
        #'     it will be recycled as necessary to match the number of columns
        #'     of X when used.
        #' @param intercept A logical vector of length one; should the linear
        #'     mean include an intercept? The default is TRUE.
        initialize = function(hypers = 0, intercept = TRUE) {
            self$hypers = hypers
            self$intercept = intercept
        }
    )
)
