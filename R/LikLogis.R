#' @title
#' Logistic Likelihood
#'
#' @description
#' The likelihood of \eqn{y_i} given \eqn{f_i} is given by the logistic
#' function with input \eqn{y_i f_i} (where the binary labels for \eqn{y}
#' are encoded by 1 and -1). It takes no hyperparameters.
#'
#' @details
#' The log likelihood is given by
#'
#' \deqn{
#'     \sum_i -\log \left(1 + \exp \left( -y_i f_i \right) \right)
#' }{%
#'     \sum_i -log ( 1 + exp( -y_i f_i ) )
#' }
#'
#' @export
LikLogis <- R6::R6Class(
    classname = "LikLogis",
    inherit = LikelihoodFunction,
    public = list(

        ## Data members
        #' @field name A character vector of length one giving the likelihood
        #'     function's name; "logistic"
        name = "logistic",
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
            return(plogis(y * f, log = TRUE))
        },
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
            if ( is.null(hypers) ) {
                hypers <- self$hypers
            }
            if ( length(y) != length(f) ) {
                stop("y and f must have the same length.")
            }
            p <- plogis(f)
            if ( order == 1 ) {
                return( ( (y+1) / 2 ) - p )
            } else if ( order == 2 ) {
                return( -p * (1 - p) )
            } else if ( order == 3 ) {
                d2lp <- -p * (1 - p)
                return( 2 * d2lp * (0.5 - p) )
            } else {
                stop("Derivatives of the log of the logistic likelihood ",
                     "function are provided only up to order 3.")
            }
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
