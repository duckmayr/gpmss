#' @title
#' Gaussian Likelihood
#'
#' @description
#' The likelihood of \eqn{y_i} given \eqn{f_i} is normal with mean \eqn{f_i}
#' and variance \eqn{\sigma_y^2}. It takes a single hyperparameter,
#' \eqn{\log(\sigma_y)}{log(\sigma_y)}.
#'
#' @section Warning:
#' Note the hyperparameter should be given on the log scale!
#'
#' @export
LikGauss <- R6::R6Class(
    classname = "LikGauss",
    inherit = LikelihoodFunction,
    public = list(

        ## Data members
        #' @field name A character vector of length one giving the likelihood
        #'     function's name; "Gaussian"
        name = "Gaussian",
        #' @field hypers A numeric vector giving the likelihood function's
        #'     hyperparameters
        hypers = numeric(),

        ## Methods
        #' @description
        #' Compute log probability of outcomes given function values
        #' @param y The observed outcomes
        #' @param f The function values (generally a posterior mean or mode)
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     likelihood function. If NULL (the default), the hypers data
        #'     member is used.
        lp = function(y, f, hypers = NULL) {
            if ( is.null(hypers) ) {
                hypers <- self$hypers
            }
            if ( length(y) != length(f) ) {
                stop("y and f must have the same length.")
            }
            return(dnorm(y, mean = f, sd = exp(hypers[1]), log = TRUE))
        },
        #' @description
        #' Compute partial derivatives of log likelihood function with respect
        #' to f
        #' @section Note:
        #' The \code{f_derivative()} method is currently unimplemented for
        #' the Gaussian likelihood
        #' (it is not needed for parameter optimization).
        #' @param y The observed outcomes
        #' @param f The function values (generally a posterior mean or mode)
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     likelihood function. If NULL (the default), the hypers data
        #'     member is used.
        #' @param order An integer vector of length one giving order of
        #'     derivative desired; the default is 1
        f_derivative = function(y, f, hypers = NULL, order = 1) {
            return(NULL)
        },

        ## Constructor
        #' @description
        #' Create a new LikelihoodFunction object
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     likelihood function; a numeric vector of length one giving
        #'     the (log of) the standard deviation of the normal distribution
        #'     for the likelihood
        initialize = function(hypers = numeric()) {
            self$hypers = hypers
        }
    )
)
