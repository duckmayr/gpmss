#' @title
#' Gaussian Process Regression Model for Exact Inference
#'
#' @description
#' Exact inference for GP regression; this inference method can only be used
#' with a Gaussian likelihood (see \code{\link{LikGauss}})
#'
#' @details
#' In Gaussian process (GP) regression, we consider that our outcomes \eqn{y}
#' are noisy observations of an unknown function of the predictors \eqn{X}:
#'
#' \deqn{
#'     y = f \left( X \right) + \varepsilon.
#' }{%
#'     y = f (X) + \epsilon.
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
#'     f (X) ~ N ( m (X), K (X, X) ),
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
#' For GP regression, the posterior distribution has an analytical solution:
#'
#' \deqn{
#'     f \mid y, X
#'     \sim
#'     \mathcal{N} \left(
#'         \mu + K K_y^{-1} \left( y - \mu \right),
#'         K - K K_y^{-1} K
#'     \right),
#' }{%
#'     f | y, X ~ N(\mu + K K_y^{-1} (y - \mu), K - K K_y^{-1} K),
#' }
#'
#' where \eqn{K} denotes the prior covariance function
#' (and we have suprressed in the notation that here it is evaluated at
#' \eqn{X}), \eqn{\mu} is the prior mean of \eqn{f}, and
#' \eqn{K_y^{-1} = (K + \sigma_y^2 I)^{-1}}.
#' For new test cases \eqn{X_\ast}{X*}, \eqn{f_\ast}{f*} is distributed
#'
#' \deqn{
#'     f_\ast \mid y, X, X_\ast
#'     \sim
#'     \mathcal{N} \left(
#'         \mu_\ast + K(X_\ast, X) K_y^{-1} \left( y - \mu \right),
#'         K(X_\ast) - K(X_\ast, X) K_y^{-1} K(X, X_\ast)
#'     \right),
#' }{%
#'     f* | y, X, X* ~ N(\mu* + K(X*, X) K_y^{-1} (y - \mu),
#'                      K(X*) - K(X*, X) K_y^{-1} K(X, X*)),
#' }
#'
#' where \eqn{\mu_\ast}{\mu*} is the prior mean function evaluated at
#' \eqn{X_\ast}{X*}.
#'
#' @examples
#' ## Simulate data
#' set.seed(123)
#' n <- 100
#' x <- rnorm(n, sd = 2)
#' y <- sin(x) + rnorm(n, sd = sqrt(0.1))
#' ## Train model
#' m <- GPR$new(y ~ x)
#' ## Visualize posterior over f at training points
#' plot(x, y, pch = "+", col = "#33333387")
#' curve(sin(x), lty = 2, add = TRUE)
#' o <- order(x)
#' lines(x[o], m$post_mean[o])
#' s <- sqrt(diag(m$post_cov))[o]
#' h <- m$post_mean[o] + 1.96 * s
#' l <- m$post_mean[o] - 1.96 * s
#' polygon(x = c(x[o], rev(x[o])), y = c(h, rev(l)),
#'         border = NA, col = "#87878733")
#' legend(x = -3, y = max(y) + 0.1, bty = "n",
#'        pch = c(3, NA, NA, 15), lty = c(NA, 2, 1, NA),
#'        col = c("#333333", "black", "black", "#87878733"),
#'        legend = c("Observations", "True f(x)", "Post. mean", "95% CI"))
#' ## Get negative log marginal likelihood of the model
#' m$nlml()
#' ## Look at derivatives of negative log marginal likelihood wrt hypers
#' m$dnlml()
#'
#' @export
GPR <- R6::R6Class(
    classname = "GPR",
    inherit = GPModel,
    public = list(
        ## Data members
        #' @field y The outcomes; should be a numeric vector.
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
        #' @field L A numeric matrix such that
        #'     \eqn{L L^T = K + \sigma_y^2 I}
        L = matrix(),
        #' @field alpha A numeric vector equal to
        #'     \eqn{(K + \sigma_y^2 I)^{-1} (y - \mu)}
        alpha = numeric(),
        #' @field post_mean A numeric vector giving the posterior mean at the
        #'     training data
        post_mean = numeric(),
        #' @field post_cov A numeric vector giving the posterior covariance at
        #'     the training data
        post_cov = matrix(),
        #' @field prior_mean A numeric vector giving the prior mean at the
        #'     training data
        prior_mean = numeric(),
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
        #'     This is restricted to be \code{\link{LikGauss}}.
        #'     If a generator is provided rather than an object,
        #'     an object will be created with the default value for the
        #'     hyperparameters.
        #' @param meanfun An object inheriting from class
        #'     \code{\link{MeanFunction}}, or a generator for such a
        #'     class. This is the mean function for the GP prior (see Details).
        #'     The default is \code{\link{MeanZero}}.
        #'     If a generator is provided rather than an object,
        #'     an object will be created with the default value for the
        #'     hyperparameters.
        #' @param covfun An object inheriting from class
        #'     \code{\link{CovarianceFunction}}, or a generator for such a
        #'     class.
        #'     This is the covariance function for the GP prior (see Details).
        #'     The default is \code{\link{CovSEard}}.
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
        #' @param ... Other arguments specific to a particular model; unused.
        initialize = function(formula, data,          ## data specification
                              likfun = LikGauss,      ## likelihood function
                              meanfun = MeanZero,     ## prior mean function
                              covfun = CovSEard,      ## prior cov  function
                              optimize = FALSE,       ## auto optimize hypers?
                              force_all_levs = FALSE, ## keep ALL factor levels?
                              ...) {
            ## Ensure we're working with Gaussian likelihood
            if ( inherits(likfun, "R6ClassGenerator") ) {
                likfun <- likfun$new()
            }
            if ( !inherits(likfun, "LikGauss") ) {
                stop("GPR *must* use Gaussian likelihood.")
            }
            ## Use superclass initializer, making sure to *not* pass "data"
            ## if it is not present
            if ( missing(data) ) {
                super$initialize(formula, likfun = likfun,
                                 meanfun = meanfun, covfun = covfun,
                                 optimize = optimize,
                                 force_all_levs = force_all_levs)
            } else {
                super$initialize(formula, data, likfun, meanfun, covfun,
                                 optimize, force_all_levs)
            }
            self$train()
        },
        ## Methods
        #' Train the GP model, providing a characterization of the posterior
        #' of the function of interest at the input values given in the
        #' training data.
        #' @param ... Additional arguments affecting the inference calculations
        train = function(...) {
            K <- self$covfun$cov(self$X)
            s <- exp(2 * self$likfun$hypers[1])
            L <- t(chol(K + diag(s, length(self$y))))
            m <- self$meanfun$mean(self$X)
            a <- L %//% (self$y - m)
            v <- solve(L, K)
            V <- K - crossprod(v)
            self$L <- L
            self$alpha <- a
            self$post_mean <- m + K %*% a
            self$post_cov  <- V
            self$prior_mean <- m
        },
        #' @description
        #' Characterize the posterior predictive distribution of the function
        #' of interest at new test points.
        #' @param newdata A data frame containing the data for the new test
        #'     points
        #' @param ... Additional arguments affecting the predictions produced
        predict = function(newdata, ...) {
            X     <- self$X
            Xstar <- construct_Xstar(self, newdata)
            Kss   <- self$covfun$cov(Xstar)
            Ks    <- self$covfun$cov(Xstar, X)
            fmean <- self$meanfun$mean(Xstar) + Ks %*% self$alpha
            fcov  <- solve(self$L, t(Ks))
            fcov  <- Kss - crossprod(fcov)
            sy2   <- exp(2 * self$likfun$hypers[1])
            ycov  <- fcov + diag(sy2, nrow = nrow(fcov))
            return(
                list(
                    fmean = fmean,
                    fcov  = fcov,
                    ymean = fmean,
                    ycov  = ycov
                )
            )
        },
        #' @description
        #' Caclulate the negative log marginal likelihood of the GP model.
        #' @param ... Additional arguments affecting the calculation
        nlml = function(...) {
            return( c(t(self$y - self$prior_mean) %*% self$alpha / 2
                      + sum(log(diag(self$L)))
                      + length(self$y) * log(2 * pi) / 2) )
        },
        #' @description
        #' Caclulate the gradient of the negative log marginal likelihood of
        #' the GP model with respect to the hyperparameters of the mean,
        #' covariance, and likelihood functions.
        #' @param ... Additional arguments affecting the calculation
        dnlml = function(...) {
            K  <- self$covfun$cov(self$X)
            Q <- self$L %//% diag(length(self$y)) - tcrossprod(self$alpha)
            res <- list()
            res[["mean"]] <- numeric(length(self$meanfun$hypers))
            for ( i in seq_along(res[["mean"]]) ) {
                d <- self$meanfun$parameter_derivative(self$X, param = i)
                res[["mean"]][i] <- -crossprod(d, self$alpha)
            }
            res[["cov"]]  <- numeric(length(self$covfun$hypers))
            for ( i in seq_along(res[["cov"]]) ) {
                d <- self$covfun$parameter_derivative(self$X, param = i, K = K)
                res[["cov"]][i] <- 0.5 * sum(Q * d)
            }
            res[["lik"]]  <- exp(2 * self$likfun$hypers[1]) * sum(diag(Q))
            return(res)
        },
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
        #' @param ... Other arguments affecting the calculations (unused)
        margins = function(variables = NULL, base_categories = NULL,
                           differences = NULL, indices = NULL, ci = 0.95,
                           force = FALSE, ...) {
            ## We make sure we have variables specified
            if ( is.null(variables) ) {
                variables <- 1:ncol(self$X)
            }
            ## and that they are stored as variable names rather than numbers
            if ( is.numeric(variables) ) {
                variables <- colnames(self$X)[variables]
            }
            ## If force = FALSE and we've already calculated them,
            ## we remove them from the variable list
            if ( !force ) {
                variables <- setdiff(variables, names(self$marginal_effects))
            }
            ## Separate the factor, binary, and continuous variables
            ## (fnames, bnames, and cnames respectively;
            ##  vnames is name of all vars)
            vnames  <- colnames(self$X)
            factors <- lapply(attr(self$X, "contrasts"), "rownames")
            fnames  <- unlist(sapply(seq_along(factors), function(i) {
                paste0(names(factors)[i], factors[[i]])
            }))
            bnames  <- vnames[apply(self$X, 2, function(x) {
                return(length(unique(x)) == 2)
            })]
            bnames  <- setdiff(bnames, fnames)
            bnames  <- union(bnames, names(differences))
            cnames  <- setdiff(vnames, c(fnames, bnames))
            ## Length variables
            n  <- nrow(self$X)
            m  <- length(variables)
            ## It will be useful to have the kernel
            K  <- self$covfun$cov(self$X)
            ## Get effects for each variable
            for ( d in variables ) {
                x   <- self$X[ , d]
                if ( d %in% cnames ) { ## if this is a continuous variable
                    ## Do partial derivatives of prior mean and covariance
                    i   <- which(colnames(self$X) == d)
                    md  <- self$meanfun$input_derivative(self$X, dimension = i)
                    Kd  <- self$covfun$input_derivative(self$X, dimension = i,
                                                        K = K)
                    Kdd <- self$covfun$input_derivative(self$X, dimension = i,
                                                        order = 2, K = K)
                    ## Get mean and variance of derivative of f wrt d
                    E   <- c(md + Kd %*% self$alpha)
                    v   <- solve(self$L, Kd)
                    V   <- Kdd - crossprod(v)
                    rownames(V) <- colnames(V) <- rownames(self$X)
                } else if ( d %in% bnames ) { ## if this is a binary variable
                    ## Get mean and variance of f(1) - f(0)
                    if ( d %in% names(differences) ) {
                        vals <- differences[[d]]
                    } else {
                        vals <- sort(unique(x))
                    }
                    X1      <- self$X
                    X1[,d]  <- vals[2]
                    X0      <- self$X
                    X0[,d]  <- vals[1]
                    Xs      <- rbind(X1, X0)
                    Ks      <- self$covfun$cov(Xs, self$X)
                    fs      <- Ks %*% self$alpha + self$meanfun$mean(Xs)
                    Kss     <- self$covfun$cov(Xs)
                    v       <- solve(self$L, t(Ks))
                    Kp      <- Kss - crossprod(v)
                    nn      <- nrow(Xs)
                    i1      <- 1:n
                    i2      <- (n+1):nn
                    V       <- Kp[i1, i1] + Kp[i2, i1] + Kp[i1, i2] + Kp[i2, i2]
                    rownames(V) <- colnames(V) <- rownames(self$X)
                    E       <- c(fs[i1] - fs[i2])
                } else { ## if this is a categorical variable
                    ## Figure out the baseline category
                    b0 <- names(factors)[which(sapply(names(factors),grepl,d))]
                    if ( !is.null(base_categories)
                         & b0 %in% names(base_categories)
                         & isTRUE(base_categories[[b0]] %in% factors[[b0]]) ) {
                        paste0(b0, base_categories[[b0]])
                    } else {
                        b <- paste0(b0, factors[[b0]][1])
                    }
                    ## If this *is* the baseline category, skip it
                    if ( d == b ) {
                        next()
                    }
                    ## Get mean and variance of f(d) - f(base_cat)
                    ## First create copy of X where all observations are type d
                    X1      <- self$X
                    X1[,paste0(b0, factors[[b0]])] <- 0
                    X1[,d]  <- 1
                    ## Then create copy of X where all observations are type b
                    X0      <- self$X
                    X0[,paste0(b0, factors[[b0]])] <- 0
                    X0[,b]  <- 1
                    ## Now get the posterior covariance of those observations
                    Xs      <- rbind(X1, X0)
                    Ks      <- self$covfun$cov(Xs, self$X)
                    Kss     <- self$covfun$cov(Xs)
                    v       <- solve(self$L, t(Ks))
                    Kp      <- Kss - crossprod(v)
                    ## And the posterior mean
                    fs      <- Ks %*% self$alpha + self$meanfun$mean(Xs)
                    ## And then we can get the mean and covariance of the diff
                    nn      <- nrow(Xs)
                    i1      <- 1:n
                    i2      <- (n+1):nn
                    V       <- Kp[i1, i1] + Kp[i2, i1] + Kp[i1, i2] + Kp[i2, i2]
                    rownames(V) <- colnames(V) <- rownames(self$X)
                    E       <- c(fs[i1] - fs[i2])
                }
                self$marginal_effects[[d]] <- list(
                    mean = E,
                    covariance = V
                )
            }
            ## Average the marginal effects over the requested indices
            if ( is.null(indices) ) {
                indices <- 1:length(self$y)
            }
            ames <- sapply(self$marginal_effects, function(x) {
                return( mean(x[["mean"]][indices]) )
            })
            k <- 1/(length(indices)^2)
            ame_sds <- sapply(self$marginal_effects, function(x) {
                return( k * sum(x[["covariance"]][indices, indices]) )
            })
            self$average_marginal_effects <- data.frame(
                Variable = names(ames),
                Mean = ames,
                LB = qnorm((1 - ci) / 2, mean = ames, sd = ame_sds),
                UB = qnorm(1 - ((1 - ci) / 2), mean = ames, sd = ame_sds)
            )
            invisible(self)
        }
    )
)
