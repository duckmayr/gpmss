#' @title
#' Mean Function Class
#'
#' @description
#' This provides a parent class for all mean function classes for gpmss.
#'
#' "A Gaussian process is completely specified by its mean function and
#' covariance function" (Rasmussen and Williams 2006, 13). The mean
#' function gives the prior mean for function output given input values.
#'
#' We implement the following mean functions:
#' \describe{
#'     \item{MeanZero}{A zero mean; for any set of inputs, the GP's prior
#'         mean is zero.}
#'     \item{MeanLinear}{A linear mean; denoting its hyper vector by
#'         \eqn{\beta}, for an input matrix \eqn{X}, the prior mean is
#'         \eqn{X \beta}.}
#' }
#'
#' Users may also define their own mean functions for use with gpmss
#' functions. A mean function class must provide a \code{mean()} method,
#' a \code{parameter_derivative()} method, and a \code{input_derivative()}
#' method. It must have data members \code{hypers} and \code{name}.
#'
#' The \code{mean()} method is used to actually call the mean function.
#' The first argument should be \code{X}, the first set of inputs.
#' The second argument should be \code{hypers}, a numeric vector giving all the
#' hyperparameters of the mean function; its default value should be
#' \code{NULL}, in which case, the data member \code{hypers} is utilized.
#' (We allow passing different hypers for convenience in hyperparameter
#' optimization).
#'
#' The \code{parameter_derivative()} method is used to obtain the partial
#' derivatives of the mean function with respect to each of its
#' hyperparameters. (This is useful for optimization). Its arguments should be
#' the same as the \code{mean()} method, with one additional argument,
#' \code{param}, giving the index of \code{hypers} of the parameter with
#' respect to which the derivative should be taken.
#'
#' The \code{input_derivative()} method is used to obtain the partial
#' derivatives of the mean function with respect to its inputs.
#' (This is useful for obtaining marginal effects). Its arguments should be
#' the same as the \code{mean()} method, with an additional argument:
#' \code{dimension}, an integer vector of length one giving the dimension of X
#' with respect to which the derivative is being taken.
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
MeanFunction <- R6::R6Class(
    classname = "MeanFunction",
    public = list(

        ## Data members
        #' @field name A character vector of length one giving the mean
        #'     function's name
        name = character(),
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
        mean = function(X, hypers = NULL) return(NULL),
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
            return(NULL)
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
            return(NULL)
        },

        ## Constructor
        #' @description
        #' Create a new MeanFunction object
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     mean function.
        initialize = function(hypers = numeric()) {
            self$hypers = hypers
        }
    )
)
