#' @title
#' Gaussian Process Model
#'
#' @description
#' GPModel provides a parent class for all Gaussian process (GP) model classes
#' in \code{gpmss}. A GP model class keeps track of the data and various
#' posterior quantities, as well as the inference method and likelihood,
#' mean, and covariance functions. It has a number of methods to produce,
#' summarize, or plot quantities of interest such as predictive distributions
#' or avarage marginal effects.
#'
#' @details
#' In Gaussian process (GP) regression, we consider that our outcomes \eqn{y}
#' are noisy observations of an unknown function of the predictors \eqn{X}:
#'
#' \deqn{
#'     y = f \left( X \right) + \varepsilon.
#' }{%
#'     y = f ( X ) + \varepsilon.
#' }
#'
#' We place a GP prior over the values of \eqn{f \left( X \right)}{f(X)},
#' which says that before observing \eqn{y}, we believe the values of
#' \eqn{f \left( X \right)}{f(X)} are distributed
#'
#' \deqn{
#'     f \left( X \right)
#'     \sim
#'     \mathcal{N} \left(
#'         m \left( X \right),
#'         K \left( X, X \right)
#'     \right),
#' }{%
#'     f ( X ) \sim N ( m ( X ), K ( X, X ) ),
#' }
#'
#' where \eqn{m \left( X \right)}{m ( X )} is a [mean function][MeanFunction]
#' and \eqn{K \left( X, X \right)}{K ( X, X )} is a
#' [covariance function][CovarianceFunction].
#' In other words, the prior mean over function outputs is a function of the
#' predictors, as is the prior covariance over function outputs.
#' Application of Bayes' rule in conjunction with Gaussian identities gives
#' us the posterior distribution of the values of \eqn{f\left(X\right)}{f(X)}
#' after observing \eqn{y}. This is a flexible Bayesian non-parametric
#' framework to reason about the relationship between predictors and outcomes
#' without making strong functional form assumptions.
#'
#' GP classification works similarly, though now \eqn{f\left(X\right)}{f(X)}
#' is a latent function that is pushed through a link function to produce
#' outcomes, not unlike generalized linear models. While the exact posterior
#' for GP regression can be derived analytically, in the classification case
#' the posterior is analytically untractable, so several approximations to
#' the posterior are available, and the posterior distribution can be
#' simulated via MCMC methods.
#'
#' A GPModel child class implements one such inference method for GP
#' regression or classification. The GP models already implemented by
#' \code{gpmss} are
#'
#' \describe{
#'     \item{GPR}{
#'         Exact inference for GP regression.
#'         This is only possible with a [Gaussian likelihood][LikGauss].
#'     }
#'     \item{GPCLA}{
#'         The Laplace approximation to the posterior for GP classification.
#'     }
#' }
#'
#' Users may also define their own GP models that utilize gpmss functions.
#' A GPModel subclass should provide a \code{train()} method,
#' a \code{predict()} method, a \code{nlml()} method, a \code{dnlml()}
#' method, and a \code{margins()} method.
#' A subclass may wish to provide a custome \code{optimize()} method for
#' optimizing hyperparameters, though the parent class' \code{optimize()}
#' method should be general purpose enough to work for any subclass designed
#' according to the description provided here.
#' It should have data members \code{name}, \code{y}, \code{X},
#' \code{meanfun} (to hold an instance of a subclass of [MeanFunction]),
#' \code{covfun} (to hold an instance of a subclass of [CovarianceFunction]),
#' \code{likfun} (to hold an instance of a subclass of [LikelihoodFunction]),
#' \code{marginal_effects}, and \code{average_marginal_effects},
#' but which other data members are appropriate may vary by inference method.
#' For most models, a data member \code{L} should be present,
#' which is a (lower) Cholesky decomposition of a matrix of central importance
#' in posterior inference, as well as \code{alpha}, which is usually
#' \eqn{K^{-1} \left(\mu_{post} - \mu_{prior}\right)}{%
#'      K^{-1} (\mu_{post} - \mu_{prior})},
#' where \eqn{\mu} is the posterior or prior mean of the function of interest
#' respectively.
#' But, again, one should consult specific models' help files
#' as the data members that are useful will differ between models.
#'
#' The \code{train()} method is used to characterize the posterior distribution
#' of the function of interest at the training cases provided in the
#' \code{data} argument of the constructor
#' (or found in the calling environment given the formula provided).
#' For many models, arguments will not be needed; the quantities necessary for
#' inference should be provided in the constructor. However, an ellipses
#' (\code{...}) argument should always conclude the argument list.
#'
#' The \code{predict()} method is used to characterize the posterior
#' predictive distribution of the function of interest at new test cases.
#' The first argument should be \code{newdata},
#' a data frame containing the data for the new test cases.
#' For many models, these may be the only arguments needed,
#' but additional arguments may be appropriate.
#' An ellipses (\code{...}) argument should always conclude the argument list.
#'
#' The \code{nlml()} method is used to calculate the negative log marginal
#' likelihood.
#' For many models, arguments will not be needed; the quantities necessary for
#' inference should be provided in the constructor or generated during the
#' call to \code{train()}. However, an ellipses
#' (\code{...}) argument should always conclude the argument list.
#'
#' The \code{dnlml()} method is used to calculate the gradient of the negative
#' log marginal likelihood with respect to the hyperparameters of the mean,
#' covariance, and likelihood functions.
#' (This is used for hyperparameter optimization).
#' For many models, arguments will not be needed; the quantities necessary for
#' inference should be provided in the constructor or generated during the
#' call to \code{train()}. However, an ellipses
#' (\code{...}) argument should always conclude the argument list.
#'
#' The \code{margins()} method is used to calculate marginal (and/or partial)
#' effects of predictors.
#' The arguments may be method dependent.
#' An ellipses (\code{...}) argument should always conclude the argument list.
#'
#' The data member \code{name} should be hard-coded within the class definition;
#' it is used for printing purposes (potentially in other functions). It should
#' be a public member so that it can be accessed directly.
#'
#' @export
GPModel <- R6::R6Class(
    classname = "GPModel",
    public = list(
        ## Data members
        #' @field y The outcomes; should be a numeric vector.
        #'     Depending on the model there may be further restrictions
        #'     on the nature of the data.
        #'     This field is usually generated automatically during object
        #'     construction and does not generally need to be interacted
        #'     with directly by the user.
        y = numeric(),
        #' @field X The predictors; should be a numeric matrix.
        #'     This field is usually generated automatically during object
        #'     construction and does not generally need to be interacted
        #'     with directly by the user.
        X = matrix(),
        #' @field terms The \code{\link[stats]{terms}} object related to the
        #'     formula provided in model construction (this is useful for
        #'     later prediction).
        terms = NULL,
        #' @field force_all_levs A logical vector of length one recording the
        #'     user's preference for dropping unused factor levels
        force_all_levs = FALSE,
        #' @field meanfun The prior mean function; must be an object of class
        #'     \code{\link{MeanFunction}}
        meanfun   = NULL,
        #' @field covfun The prior covariance function; must be an object of
        #'     class \code{\link{CovarianceFunction}}
        covfun    = NULL,
        #' @field likfun The likelihood function; must be an object of class
        #'     \code{\link{LikelihoodFunction}}
        likfun    = NULL,
        #' @field marginal_effects A list with an element for each predictor in
        #'     the model marginal effects have been requested for (via the
        #'     \code{margins()} method). Each element is a list of length two,
        #'     with an element "mean", each entry i of which gives the mean of
        #'     the distribution of the marginal effect of the predictor on the
        #'     ith observation, and an element "covariance", giving the
        #'     covariance matrix for the distribution of marginal effects of
        #'     that predictor.
        marginal_effects = list(),
        #' @field average_marginal_effects A dataframe with a row for each
        #'     predictor in the model that marginal effects have been requested
        #'     for (via the \code{margins()} method) and columns:
        #'     "Variable", giving the name of the predictor;
        #'     "Mean", giving the mean of the average marginal effect of that
        #'     predictor;
        #'     "LB", giving the lower bound on the requested confidence
        #'     interval (see the \code{margins()} method);
        #'     and "UB", giving the upper bound on the CI.
        average_marginal_effects = data.frame(),

        ## Constructor
        #' @description
        #' Create a new GPModel object
        #' @param formula A formula object giving the variable name of the
        #'     outcomes on the left-hand side and the predictors on the
        #'     right-hand side, a la \code{\link[stats]{lm}}
        #' @param data An optional data frame where the variables in
        #'     \code{formula} are to be found. If not found there,
        #'     we search for the variables elsewhere, generally the
        #'     calling environment.
        #' @param likfun An object inheriting from class
        #'     \code{\link{LikelihoodFunction}}, or a generator for such a
        #'     class. This is the likelihood function for the GP model.
        #'     If a generator is provided rather than an object,
        #'     an object will be created with the default value for the
        #'     hyperparameters.
        #' @param meanfun An object inheriting from class
        #'     \code{\link{MeanFunction}}, or a generator for such a
        #'     class. This is the mean function for the GP prior (see Details).
        #'     If a generator is provided rather than an object,
        #'     an object will be created with the default value for the
        #'     hyperparameters.
        #' @param covfun An object inheriting from class
        #'     \code{\link{CovarianceFunction}}, or a generator for such a
        #'     class.
        #'     This is the covariance function for the GP prior (see Details).
        #'     If a generator is provided rather than an object,
        #'     an object will be created with the default value for the
        #'     hyperparameters.
        #' @param optimize A logical vector of length one; if \code{TRUE},
        #'     the hyperparameters of the mean, covariance, and likelihood
        #'     functions are optimized automatically as part of the model
        #'     construction process (see Details).
        #'     The default is \code{FALSE}, meaning the hyperparameters given
        #'     at the time of creating those objects are used for inference.
        #' @param force_all_levs A logical vector of length one; if \code{TRUE},
        #'     unused factor levels in right-hand side variables are
        #'     \strong{not} dropped. The default is \code{FALSE}.
        #' @param ... Other arguments specific to a particular model;
        #'     see the help file for specific models for more details.
        initialize = function(formula, data,          ## data specification
                              likfun = NULL,          ## likelihood function
                              meanfun = NULL,         ## prior mean function
                              covfun = NULL,          ## prior cov  function
                              optimize = FALSE,       ## auto optimize hypers?
                              force_all_levs = FALSE, ## keep ALL factor levels?
                              ...) {
            ## Fix mean, cov, & lik funcs
            if ( inherits(meanfun, "R6ClassGenerator") ) {
                meanfun <- meanfun$new()
            }
            if ( inherits(meanfun, "MeanFunction") ) {
                self$meanfun <- meanfun
            } else {
                stop("Provide an object inheriting from MeanFunction, ",
                     "or an R6ClassGenerator for one, for meanfun.")
            }
            if ( inherits(covfun, "R6ClassGenerator") ) {
                covfun <- covfun$new()
            }
            if ( inherits(covfun, "CovarianceFunction") ) {
                self$covfun <- covfun
            } else {
                stop("Provide an object inheriting from CovarianceFunction, ",
                     "or an R6ClassGenerator for one, for covfun.")
            }
            if ( inherits(likfun, "R6ClassGenerator") ) {
                likfun <- likfun$new()
            }
            if ( inherits(likfun, "LikelihoodFunction") ) {
                self$likfun <- likfun
            } else {
                stop("Provide an object inheriting from LikelihoodFunction, ",
                     "or an R6ClassGenerator for one, for likfun.")
            }
            ## Prepare the data
            ## Fix the formula for one-hot encoding
            formula <- convert_formula(formula)
            ## The next few lines were adapted from stats::lm,
            ## with some additions to deal with variables that should be
            ## treated like factors but are stored as character vectors
            cl <- match.call()
            mf <- cl[c(1L, match(c("formula", "data"),names(cl),0L))]
            mf[[1L]] <- quote(stats::model.frame)
            names(mf)[2] <- "formula"
            if ( force_all_levs ) {
                mf$drop.unused.levels <- FALSE
            } else {
                mf$drop.unused.levels <- TRUE
            }
            self$force_all_levs <- force_all_levs
            mf <- eval(mf, parent.frame())
            ## It will be useful in prediction to store the terms
            self$terms <- attr(mf, "terms")
            ## The 1st column of mf will be the outcomes, the others the inputs
            self$y <- c(as.matrix(mf[ , 1]))
            ## Get some info to deal with factors
            ff <- sapply(mf, function(x) is.factor(x) | is.character(x))
            cc <- list()
            for ( i in seq_along(ff) ) {
                if ( ff[i] ) {
                    this_levs  <- levels(mf[[names(ff)[i]]])
                    if ( is.null(this_levs) ) {
                        this_levs <- unique(mf[[names(ff)[i]]])
                    }
                    this_contr <- diag(length(this_levs))
                    dimnames(this_contr) <- list(this_levs, this_levs)
                    cc[[length(cc)+1]] <- this_contr
                }
            }
            names(cc) <- names(ff)[ff]
            self$X <- stats::model.matrix(formula, mf, contrasts.arg = cc)
            if ( optimize ) {
                warning("Hyper optimization not yet implemented.")
            }
        },
        ## Methods
        #' Train the GP model, providing a characterization of the posterior
        #' of the function of interest at the input values given in the
        #' training data.
        #' @param ... Additional arguments affecting the inference calculations
        train = function(...) return(NULL),
        #' @description
        #' Characterize the posterior predictive distribution of the function
        #' of interest at new test points.
        #' @param newdata A data frame containing the data for the new test
        #'     points
        #' @param ... Additional arguments affecting the predictions produced
        predict = function(newdata, ...) return(NULL),
        #' @description
        #' Caclulate the negative log marginal likelihood of the GP model.
        #' @param ... Additional arguments affecting the calculation
        nlml = function(...) return(NULL),
        #' @description
        #' Caclulate the gradient of the negative log marginal likelihood of
        #' the GP model with respect to the hyperparameters of the mean,
        #' covariance, and likelihood functions.
        #' @param ... Additional arguments affecting the calculation
        dnlml = function(...) return(NULL),
        #' @description
        #' Caclulate marginal effects of predictors.
        #' @details
        #' The derivative of a GP is also a GP, giving us easy access to the
        #' distribution of marginal effects of predictors. The first time this
        #' method is called, it calculates and stores the distribution of
        #' pointwise partial derivatives of \eqn{f} for each specified predictor
        #' in the model (if no predictors are specified, marginal effects are
        #' calculated for all predictors). If a predictor is binary, instead
        #' the distribution of the difference in \eqn{f} between the predictor
        #' taking the binary variable's highest value and its lowest value is
        #' calculated. If a predictor is categorical, a similar calculation is
        #' made, but comparing each of the labels except one to a baseline
        #' label. The user may specify to use a similar calculation for two
        #' user specified values for continuous variables rather than using
        #' the partial derivatives (in some cases this may be more easily
        #' interpretable). Additional calls to \code{margins()} may calculate
        #' marginal effects for predictors whose marginal effects had not yet
        #' been requested, but marginal effects for variables that have already
        #' been calculated are not re-calculated unless \code{force = TRUE}.
        #'
        #' Every time the method is called, it stores a dataframe of average
        #' marginal effects; it calculates the mean and variance of the
        #' marginal effect of each predictor over all observations, or over
        #' specified indices of observations if provided.
        #' @param variables A vector specifying the variables marginal effects
        #'     are desired for; can be an integer vector, giving the column
        #'     indices of X to get effects for, or a character vector;
        #'     if NULL (the default), effects are derived for all variables.
        #' @param base_categories If X contains contrasts,
        #'     the marginal effects will be the differences between the levels
        #'     and a base level. By default, the base level is the lowest
        #'     factor level of the contrast, but you can pass a named list to
        #'     change the base level for some or all of the variables assigned
        #'     to a contrast scheme.
        #' @param differences A named list of 2-length numeric vectors
        #'     may be provided giving the low (first element of each vector)
        #'     and high (second element) values to calculate the effect at for
        #'     continuous variables. Any elements giving values for binary or
        #'     categorical variables are ignored (this is meaningless for
        #'     binary variables as there are only two values to choose from,
        #'     and categorical variables should be controlled with the option
        #'     base_categories).
        #'     If NULL (the default), derivatives are used for marginal
        #'     effects of continuous predictors.
        #' @param indices A numeric vector of indices over which to average
        #'     marginal effects. If NULL (the default), all observations are
        #'     used.
        #' @param ci A numeric vector of length one such that 0 < ci < 1
        #'     giving the width of the confidence interval for the average
        #'     marginal effects. The default is 0.95, corresponding to a
        #'     95% confidence interval.
        #' @param force A logical vector of length one; should marginal effects
        #'     for variables who have already had effects calculated be
        #'     re-calculated? The default is FALSE. (This is useful in case the
        #'     user has requested effects for continuous variables be
        #'     calculated as differences but would now prefer derivatives,
        #'     or vice versa).
        #' @param ... Other arguments affecting the calculations
        margins = function(variables = NULL, base_categories = NULL,
                           differences = NULL, indices = NULL, ci = 0.95,
                           force = FALSE, ...) {
            return(NULL)
        },
        #' @description
        #' Set hyperparameters by optimizing the log marginal likelihood
        #' @details
        #' The default settings of this function should result in an
        #' optimization routine that is very similar to that employed in the
        #' MATLAB/Octave software GPML. For details on the optimization options,
        #' see \code{\link[mize]{mize}}
        #' @param method Optimization method; default is "CG"
        #' @param line_search Type of line search; default is "Rasmussen"
        #' @param step0 Line search value for the first step; default is
        #'     "rasmussen"
        #' @param ... Other arguments passed to \code{\link[mize]{mize}}
        optimize = function(method = "CG", line_search = "Rasmussen",
                            step0 = "rasmussen", ...) {
            funlist <- list(
                fn = self$nlml,
                gr = function(par) {
                    n1 <- length(self$meanfun$hypers)
                    n2 <- length(self$covfun$hypers)
                    n3 <- length(self$likfun$hypers)
                    if ( n1 > 0 ) {
                        self$meanfun$hypers <- par[1:n1]
                    }
                    if ( n2 > 0 ) {
                        self$covfun$hypers  <- par[(n1+1):(n1+n2)]
                    }
                    if ( n3 > 0 ) {
                        self$likfun$hypers  <- par[(n1+n2+1):(n1+n2+n3)]
                    }
                    self$train()
                    g  <- self$dnlml()
                    return(c(g$mean, g$cov, g$lik))
                },
                fg = function(par) {
                    n1 <- length(self$meanfun$hypers)
                    n2 <- length(self$covfun$hypers)
                    n3 <- length(self$likfun$hypers)
                    if ( n1 > 0 ) {
                        self$meanfun$hypers <- par[1:n1]
                    }
                    if ( n2 > 0 ) {
                        self$covfun$hypers  <- par[(n1+1):(n1+n2)]
                    }
                    if ( n3 > 0 ) {
                        self$likfun$hypers  <- par[(n1+n2+1):(n1+n2+n3)]
                    }
                    self$train()
                    g  <- self$dnlml()
                    return(
                        list(
                            fn = self$nlml(),
                            gr = c(g$mean, g$cov, g$lik)
                        )
                    )
                }
            )
            h  <- c(self$meanfun$hypers, self$covfun$hypers, self$likfun$hypers)
            h  <- mize::mize(h, funlist, method = method, step0 = step0,
                             line_search = line_search, ...)
            if ( grepl(pattern = "max", x = h$terminate$what) ) {
                warning("Convergence was not reached in optimization; ",
                        "additional calls to optimize() and/or changes to ",
                        "optimization options may be desirable.")
            }
            n1 <- length(self$meanfun$hypers)
            n2 <- length(self$covfun$hypers)
            n3 <- length(self$likfun$hypers)
            if ( n1 > 0 ) {
                self$meanfun$hypers <- h$par[1:n1]
            }
            if ( n2 > 0 ) {
                self$covfun$hypers  <- h$par[(n1+1):(n1+n2)]
            }
            if ( n3 > 0 ) {
                self$likfun$hypers  <- h$par[(n1+n2+1):(n1+n2+n3)]
            }
            invisible(self)
        }
    )
)
