#' @title
#' Zero Mean Function
#'
#' @description
#' A zero mean; for any set of inputs, this function specifies that the
#' Gaussian process's mean is zero.
#'
#' @details
#' As this function specifies a constant prior mean, it has no hyperparameters.
#'
#' @export
MeanZero <- R6::R6Class(
    classname = "MeanZero",
    inherit = MeanFunction,
    public = list(

        ## Data members
        #' @field name A character vector of length one giving the mean
        #'     function's name; "zero"
        name = "zero",
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
        mean = function(X, hypers = NULL) return(rep(0, nrow(X))),
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
            return(numeric())
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
            return(rep(0, nrow(X)))
        },

        ## Constructor
        #' @description
        #' Create a new MeanZero object
        #' @param hypers A(n empty) numeric vector giving hyperparameters for
        #'     the mean function.
        initialize = function(hypers = numeric()) {
            self$hypers = hypers
        }
    )
)
