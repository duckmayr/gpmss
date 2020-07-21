#' @title
#' Covariance Function Class
#'
#' @description
#' This provides a parent class for all covariance function classes for gpmss.
#'
#' "A Gaussian process is completely specified by its mean function and
#' covariance function" (Rasmussen and Williams 2006, 13). The covariance
#' function gives the prior covariance between function outputs given two
#' sets of input values. A valid covariance function must be a positive
#' semi-definite kernel; see Rasmussen and Williams (2006), Chapter 4 for
#' details.
#'
#' We implement the following covariance functions:
#' \describe{
#'     \item{CovSEiso}{The isometric squared exponential covariance function.
#'         This is also sometimes called the Gaussian kernel, the squared
#'         exponential kernel, or the radial basis function.}
#'     \item{CovSEard}{The squared exponential covariance function with
#'         automatic relevance determination. This is also sometimes called the
#'         anisotropic squared exponential kernel or covariance function.}
#' }
#'
#' Users may also define their own covariance functions for use with gpmss
#' functions. A covariance function class must provide a \code{cov()} method,
#' a \code{parameter_derivative()} method, and a \code{input_derivative()}
#' method. It must have data members \code{hypers} and \code{name}.
#'
#' The \code{cov()} method is used to actually call the covariance function.
#' The first argument should be \code{X}, the first set of inputs.
#' The second argument should be \code{Z}, the second set of inputs;
#' its default value should be \code{X}.
#' The third argument should be \code{hypers}, a numeric vector giving all the
#' hyperparameters of the covariance function; its default value should be
#' \code{NULL}, in which case, the data member \code{hypers} is utilized.
#' (We allow passing different hypers for convenience in hyperparameter
#' optimization).
#'
#' The \code{parameter_derivative()} method is used to obtain the partial
#' derivatives of the covariance function with respect to each of its
#' hyperparameters. (This is useful for optimization). Its arguments should be
#' the same as the \code{cov()} method, with one additional argument,
#' \code{param}, giving the index of \code{hypers} of the parameter with
#' respect to which the derivative should be taken.
#'
#' The \code{input_derivative()} method is used to obtain the partial
#' derivatives of the covariance function with respect to its inputs.
#' (This is useful for obtaining marginal effects). Its arguments should be
#' the same as the \code{cov()} method, with two additional arguments:
#' \code{dimension}, an integer vector of length one giving the dimension of X
#' with respect to which the derivative is being taken; and \code{order},
#' an integer vector of length one giving the order of partial derivative
#' desired.
#'
#' Its constructor should take one mandatory argument, \code{hypers},
#' which will be stored as the data member \code{hypers}.
#' The data member \code{hypers} should be a public member so that it can be
#' accessed and modified directly.
#'
#' The data member \code{name} should be hard-coded within the class definition;
#' it is used for printing purposes (potentially in other functions). It should
#' be a public member so that it can be accessed directly.
#'
#' @export
CovarianceFunction <- R6::R6Class(
    classname = "CovarianceFunction",
    public = list(

        ## Data members
        #' @field name A character vector of length one giving the covariance
        #'     function's name
        name = character(),
        #' @field hypers A numeric vector giving the covariance function's
        #'     hyperparameters
        hypers = numeric(),

        ## Methods
        #' @description
        #' Compute function covariance
        #' @param X The first set of input values (should be a numeric matrix)
        #' @param Z The second set of input values (should be a numeric matrix);
        #'     The default is Z = X.
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     covariance function. If NULL (the default), the hypers data
        #'     member is used.
        cov = function(X, Z = X, hypers = NULL) return(NULL),
        #' @description
        #' Compute partial derivatives of covariance function with respect to
        #' its hyperparameters
        #' @param X The first set of input values (should be a numeric matrix)
        #' @param Z The second set of input values (should be a numeric matrix);
        #'     The default is Z = X.
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     covariance function. If NULL (the default), the hypers data
        #'     member is used.
        #' @param param An integer vector of length one; which element of
        #'     \code{hypers} should the derivative be taken with respect to?
        #'     The default is 1
        parameter_derivative = function(X, Z = X, hypers = NULL, param = 1) {
            return(NULL)
        },
        #' @description
        #' Compute partial derivatives of covariance function with respect to
        #' its inputs
        #' @param X The first set of input values (should be a numeric matrix)
        #' @param Z The second set of input values (should be a numeric matrix);
        #'     The default is Z = X.
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     covariance function. If NULL (the default), the hypers data
        #'     member is used.
        #' @param dimension an integer vector of length one giving the dimension
        #'     of X with respect to which the derivative is being taken; the
        #'     default is 1
        #' @param order An integer vector of length one indicating whether the
        #'     first partial derivative (order = 1) is desired, or the cross
        #'     partial (order = 2); the default is 1
        input_derivative = function(X, Z = X, hypers = NULL,
                                    dimension = 1, order = 1) {
            return(NULL)
        },

        ## Constructor
        #' @description
        #' Create a new CovarianceFunction object
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     covariance function.
        initialize = function(hypers = numeric()) {
            self$hypers = hypers
        }
    )
)
