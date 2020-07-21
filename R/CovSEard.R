#' @title
#' Squared Exponential Covariance Function with Automatic Relevance
#' Determination
#'
#' @description
#' The anisometric squared exponential covariance function, also called the
#' squared exponential covariance function with automatic relevance
#' determination.
#'
#' @details
#' The anisometric squared exponential covariance function has D + 1
#' hyperparameters (where D is the number of dimensions in X),
#' \eqn{\sigma_f}, the scale factor, and \eqn{\boldsymbol\ell},
#' a vector of characteristic length scales.
#' The scale factor governs (on average) how far from the mean the function
#' values can be, while the length scales govern how quickly the
#' function can change given movement across a particular dimension of X;
#' or, in other words, as the function output covariance is given as a
#' function of distance in the covariate space,
#' the length scale governs what "closeness" means,
#' for each dimension.
#' This can also be interpreted as determining the relative importance of
#' dimensions, particularly for dichotomous dimensions of X
#' (hence the moniker "automatic relevance determination").
#'
#' The covariance between \eqn{f(\mathbf{x}_i)} and \eqn{f(\mathbf{x}_j)}
#' is given by
#'
#' \deqn{
#'     k \left( \mathbf{x}_i, \mathbf{x}_j \right)
#'     =
#'     \sigma_f^2
#'     \exp \left(
#'         \left( \mathbf{x}_i - \mathbf{x}_j \right)^T
#'         M \left( \mathbf{x}_i - \mathbf{x}_j \right)
#'     \right),
#' }
#'
#' where \eqn{M} is a matrix whose dth diagonal entry is \eqn{1 / \ell_d^2}.
#'
#' @section Warning:
#' Note that the hyperparameters should be stored on the log scale;
#' that is, you should supply the log of the scale factor and the log of the
#' length scale (in that order).
#'
#' @export
CovSEard <- R6::R6Class(
    classname = "CovSEard",
    inherit = CovarianceFunction,
    public = list(

        ## Data members
        #' @field name A character vector of length one giving the covariance
        #'     function's name; "anisotropic squared exponential"
        name = "anisotropic squared exponential",
        #' @field hypers A numeric vector giving the covariance function's
        #'     hyperparameters; a vector of length D+1 giving the log of the
        #'     scale factor and the log of the length scale, in that order
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
        cov = function(X, Z = X, hypers = NULL) {
            if ( is.null(hypers) ) {
                hypers <- self$hypers
            }
            if ( length(hypers) != (ncol(X)+1) ) {
                stop("CovSEiso should be called with D+1 hyperparameters.")
            }
            s <- exp(2 * hypers[1])
            M <- diag(1 / exp(2 * hypers[-1]))
            n <- nrow(X)
            m <- nrow(Z)
            K <- matrix(NA_real_, nrow = n, ncol = m)
            for ( j in 1:m ) {
                for ( i in 1:n ) {
                    d <- X[i, ] - Z[j, ]
                    K[i, j] <- s * exp(-0.5 * t(d) %*% M %*% d)
                }
            }
            return(K)
        },
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
        #'     If 1 (the default), the derivative is taken with respect to the
        #'     (log of the) scale factor; if 2 or more, it is taken with
        #'     respect to the param - 1 element of the
        #'     (log of the) length scale.
        parameter_derivative = function(X, Z = X, hypers = NULL, param = 1) {
            if ( is.null(hypers) ) {
                hypers <- self$hypers
            }
            D  <- ncol(X)
            if ( length(hypers) != (D+1) ) {
                stop("CovSEard should be called with D+1 hyperparameters.")
            }
            K  <- self$cov(X, Z, hypers)
            if ( param == 1 ) {
                return(2 * K) ## partial wrt scale factor
            } else if ( param > 1 & param < D + 2 ) {
                sqd <- outer(X[ , param-1], Z[ , param-1], "-")
                sqd <- (sqd * sqd) / exp(2 * hypers[param])
                return(sqd * K)
            } else {
                stop("param must be between 1 and D + 1")
            }
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
            if ( is.null(hypers) ) {
                hypers <- self$hypers
            }
            if ( length(hypers) != (ncol(X)+1) ) {
                stop("CovSEard should be called with D+1 hyperparameters.")
            }
            K  <- self$cov(X, Z, hypers)
            n  <- nrow(K)
            m  <- ncol(K)
            dK <- matrix(NA_real_, n, m)
            for ( j in 1:m ) {
                for ( i in 1:n ) {
                    dK[i, j] <- X[j, dimension] - X[i, dimension]
                }
            }
            ell2 <- exp(2 * hypers[dimension+1])
            if ( order == 1 ) {
                dK <- K * (dK / ell2)
            } else if ( order == 2 ) {
                dK <- K * ( (1/ell2) + ( (dK * -dK) / (ell2*ell2) ) )
            } else {
                stop("CovSEard derivatives are only programmed up to order 2")
            }
            return(dK)
        },

        ## Constructor
        #' @description
        #' Create a new CovSEiso object
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     covariance function; a vector of length two giving the log of
        #'     the scale factor and the log of the length scale, in that order
        initialize = function(hypers = numeric()) {
            self$hypers = hypers
        }
    )
)
