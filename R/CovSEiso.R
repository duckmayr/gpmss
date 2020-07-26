#' @title
#' Isotropic Squared Exponential Covariance Function
#'
#' @description
#' The isometric squared exponential covariance function, also  called the
#' Gaussian kernel, the squared exponential kernel, or the radial basis
#' function.
#'
#' @details
#' The isometric squared exponential covariance function has two
#' hyperparameters, \eqn{\sigma_f}, the scale factor, and \eqn{\ell}{ell},
#' the characteristic length scale. The scale factor governs (on average)
#' how far from the mean function values can be, while the length scale governs
#' how quickly the function can change; or, in other words, as the function
#' output covariance is given as a function of distance in the covariate space,
#' the length scale governs what "closeness" means.
#'
#' The covariance between \eqn{f(\mathbf{x}_i)}{f(x_i)}
#' and \eqn{f(\mathbf{x}_j)}{f(x_j)}
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
#' }{%
#'     k ( x_i, x_j ) = \sigma_f^2 exp [ ( x_i - x_j )^T M ( x_i - x_j ) ],
#' }
#'
#' where \eqn{M} is a matrix whose diagonal entries are
#' \eqn{1 / \ell^2}{1 / ell^2}.
#'
#' @section Warning:
#' Note that the hyperparameters should be stored on the log scale;
#' that is, you should supply the log of the scale factor and the log of the
#' length scale (in that order).
#'
#' @export
CovSEiso <- R6::R6Class(
    classname = "CovSEiso",
    inherit = CovarianceFunction,
    public = list(

        ## Data members
        #' @field name A character vector of length one giving the covariance
        #'     function's name; "isotropic squared exponential"
        name = "isotropic squared exponential",
        #' @field hypers A numeric vector giving the covariance function's
        #'     hyperparameters; a vector of length two giving the log of the
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
            if ( length(hypers) != 2 ) {
                stop("CovSEiso should be called with two hyperparameters.")
            }
            n <- nrow(X)
            m <- nrow(Z)
            s <- exp(2 * hypers[1])
            M <- diag(1 / exp(2 * hypers[2]), nrow = ncol(X))
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
        #'     (log of the) scale factor; if 2, it is taken with respect to the
        #'     (log of the) length scale.
        #' @param K An optional provision of the pre-computed kernel;
        #'     this is useful if parameter_derivative() will be called
        #'     repeatedly (for the different hypers) without the kernel
        #'     itself changing
        parameter_derivative = function(X, Z = X, hypers = NULL, param = 1,
                                        K = NULL) {
            if ( is.null(hypers) ) {
                hypers <- self$hypers
            }
            if ( length(hypers) != 2 ) {
                stop("CovSEiso should be called with two hyperparameters.")
            }
            if ( is.null(K) ) {
                K  <- self$cov(X, Z, hypers)
            }
            if ( param == 1 ) {
                return(2 * K) ## partial wrt scale factor
            } else {
                n  <- nrow(X)
                m  <- nrow(Z)
                l2 <- exp(2 * hypers[2])
                dK <- matrix(0, n, m)
                for ( j in 1:m ) {
                    for ( i in 1:n ) {
                        dK[i, j] <- sum((X[i, ] - Z[j, ])^2) / l2
                    }
                }
                return(K * dK) ## partial wrt length scale
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
        #' @param K An optional provision of the pre-computed kernel;
        #'     this is useful if parameter_derivative() will be called
        #'     repeatedly (for the different hypers) without the kernel
        #'     itself changing
        input_derivative = function(X, Z = X, hypers = NULL,
                                    dimension = 1, order = 1,
                                    K = NULL) {
            if ( is.null(hypers) ) {
                hypers <- self$hypers
            }
            if ( length(hypers) != 2 ) {
                stop("CovSEiso should be called with two hyperparameters.")
            }
            if ( is.null(K) ) {
                K  <- self$cov(X, Z, hypers)
            }
            n  <- nrow(K)
            m  <- ncol(K)
            dK <- matrix(NA_real_, n, m)
            for ( j in 1:m ) {
                for ( i in 1:n ) {
                    dK[i, j] <- X[j, dimension] - X[i, dimension]
                }
            }
            ell2 <- exp(2 * hypers[2])
            if ( order == 1 ) {
                dK <- K * (dK / ell2)
            } else if ( order == 2 ) {
                dK <- K * ( (1/ell2) + ( (dK * -dK) / (ell2*ell2) ) )
            } else {
                stop("CovSEiso derivatives are only programmed up to order 2")
            }
            return(dK)
        },

        ## Constructor
        #' @description
        #' Create a new CovSEiso object
        #' @param hypers A numeric vector giving hyperparameters for the
        #'     covariance function; a vector of length two giving the log of the
        #'     scale factor and the log of the length scale, in that order
        initialize = function(hypers = c(0, 0)) {
            self$hypers = hypers
        }
    )
)
