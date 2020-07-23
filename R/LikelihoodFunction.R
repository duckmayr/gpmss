#' @title
#' Likelihood Function Class
#'
#' @description
#' This provides a parent class for all likelihood function classes for gpmss.
#'
#' @details
#' We implement the following likelihood functions:
#' \describe{
#'     \item{LikGauss}{The likelihood of \eqn{y_i} given \eqn{f_i} is normal
#'         with mean \eqn{f_i} and variance \eqn{\sigma_y^2}.}
#'     \item{LikLogis}{The likelihood of \eqn{y_i} given \eqn{f_i} is given by
#'         the logistic function with input \eqn{y_i f_i} (where the binary
#'         labels for \eqn{y} are encoded by 1 and -1).}
#' }
#'
#' Users may also define their own likelihood functions for use with gpmss
#' functions. A likelihood function class must provide a \code{lp()} method
#' and a \code{f_derivative()} method.
#' It must have data members \code{hypers} and \code{name}.
#'
#' The \code{lp()} method calculates the log probability of observations.
#' The first argument should be \code{y}, the outcomes.
#' The second argument should be \code{f}, the function values.
#' The third argument should be \code{hypers}, a numeric vector giving all the
#' hyperparameters of the likelihood function; its default value should be
#' \code{NULL}, in which case, the data member \code{hypers} is utilized.
#' (We allow passing different hypers for convenience in hyperparameter
#' optimization).
#'
#' The \code{f_derivative()} method is used to obtain the partial
#' derivatives of the log likelihood function with respect to f.
#' (This is useful for optimization in some inference methods).
#' Its arguments should be the same as the \code{lp()} method,
#' with one additional argument, \code{order},
#' giving the order of derivative desired.
#'
#' Its constructor should take at least one argument, \code{hypers},
#' which will be stored as the data member \code{hypers}.
#' This should have a sane default, as many users may pass the class generator
#' to the constructor for a \code{\link{GPModel}} without explicitly
#' specifying hyperparameters. This should not cause much difficulty if the
#' user optimizes the hyperparameters, but a sane default should be provided
#' nonetheless.
#' The data member \code{hypers} should be a public member so that it can be
#' accessed and modified directly.
#'
#' The data member \code{name} should be hard-coded within the class definition;
#' it is used for printing purposes (potentially in other functions). It should
#' be a public member so that it can be accessed directly.
#'
#' @export
LikelihoodFunction <- R6::R6Class(
    classname = "LikelihoodFunction",
    public = list(

        ## Data members
        #' @field name A character vector of length one giving the likelihood
        #'     function's name
        name = character(),
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
        lp = function(y, f, hypers = NULL) return(NULL),
        #' @description
        #' Compute partial derivatives of log likelihood function with respect
        #' to f
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
        #'     likelihood function.
        initialize = function(hypers = numeric()) {
            self$hypers = hypers
        }
    )
)
