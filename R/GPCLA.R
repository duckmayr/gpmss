#' @title
#' Gaussian Process Classification Model Using the Laplace Approximation
#'
#' @description
#' Approximate inference for GP classification using the Laplace approximation;
#' this inference method can be used with a likelihood that accepts only values
#' -1 and 1, such as \code{\link{LikLogis}}.
#'
#' @details
#' In Gaussian process (GP) classification, we consider the likelihood of
#' outcomes \eqn{y \in \{-1, 1\}} to be
#'
#' \deqn{
#'     \prod_i \sigma \left( y_i f \left( x_i \right) \right),
#' }{%
#'     \prod_i \sigma ( y_i f(x_i) ),
#' }
#'
#' where \eqn{\sigma} is a function mapping the reals to \eqn{[0, 1]}.
#' In many dichotomous models in R, the dichotomous values are constrained to
#' be in \{0, 1\}; we follow the standard practice in GP models of instead using
#' \{-1, 1\} as this is more convenient and efficient for likelihood-related
#' computation. If \eqn{y} is provided with some other set of dichotomous
#' values, the outcomes are transformed to \{-1, 1\}, with a warning
#' (see the \code{formula} argument of the \code{new()} method for details).
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
#'
#' Unlike with [GP regression][GPR], we cannot simply apply Bayes' rule
#' and derive an anlyitical posterior. However, we can use Laplace's method
#' to approximate the posterior with a normal distribution:
#'
#' \deqn{
#'     f \mid y, X
#'     \sim
#'     \mathcal{N} \left(
#'         \mu + K \nabla \log p \left( y \mid \hat{f} \right),
#'         \left( K^{-1} + W \right)^{-1}
#'     \right),
#' }{%
#'     f | y, X ~ N(\mu + K \nabla log p(y | \hat{f}), (K^{-1} + W)^{-1} ),
#' }
#'
#' where \eqn{K} denotes the prior covariance function
#' (and we have suprressed in the notation that here it is evaluated at
#' \eqn{X}), \eqn{\mu} is the prior mean of \eqn{f}, and
#' \eqn{W = -\nabla \nabla \log p \left( y \mid \hat{f} \right)}{%
#'      W = -\nabla \nabla log p(y | \hat{f})
#' }.
#' For new test cases \eqn{X_\ast}{X*}, \eqn{f_\ast}{f*} is distributed
#'
#' \deqn{
#'     f_\ast \mid y, X, X_\ast
#'     \sim
#'     \mathcal{N} \left(
#'         \mu_\ast + K(X_\ast, X) \nabla \log p \left( y \mid \hat{f} \right),
#'         K(X_\ast) - K(X_\ast, X) (K + W^{-1})^{-1} K(X, X_\ast)
#'     \right),
#' }{%
#'     f* | y, X, X* ~ N(\mu* + \nabla log p(y | \hat{f}),
#'                      K(X*) - K(X*, X) (K + W^{-1})^{-1} K(X, X*)),
#' }
#'
#' where \eqn{\mu_\ast}{\mu*} is the prior mean function evaluated at
#' \eqn{X_\ast}{X*}.
#' (See Rasmussen and Williams (2006), Section 3.4 for details).
#'
#' @export
GPCLA <- R6::R6Class(
    classname = "GPCLA",
    inherit = GPModel,
    public = list(
        ## Data members
        #' @field y The outcomes; should be a numeric vector.
        #'     The outcome variable should only have -1 and 1 values.
        #'     If the outcome has only two values that are not 1 and -1,
        #'     the highest value is reassigned the value 1 and the lowest
        #'     value is reassigned the value -1, with a warning.
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
        #'     \eqn{L L^T = I + W^{1/2} K W^{1/2}}
        L = matrix(),
        #' @field alpha A numeric vector equal to
        #'     \eqn{\nabla \log p \left( y \mid \hat{f} \right)}{%
        #'          \nabla log p(y | \hat{f})
        #'     }
        alpha = numeric(),
        #' @field sW A numeric matrix equal to \eqn{W^{1/2}}
        sW = matrix(),
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
        #'     \code{margins()} method). Each element is itself a list.
        #'     If marginal effects on the latent function scale have been
        #'     requested, this is a list of length two,
        #'     with an element "mean", each entry i of which gives the mean of
        #'     the distribution of the marginal effect of the predictor on the
        #'     ith observation, and an element "covariance", giving the
        #'     covariance matrix for the distribution of marginal effects of
        #'     that predictor.
        #'     If marginal effects on the probability scale have been requested,
        #'     this is a list of length one with element "draws", giving the
        #'     draws of the requested marginal effects (see the Details
        #'     section of the margins() method).
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
        #'     right-hand side, a la \code{\link[stats]{lm}}.
        #'     Note the outcome variable should only have -1 and 1 values.
        #'     If the outcome has only two values that are not 1 and -1,
        #'     the highest value is reassigned the value 1 and the lowest
        #'     value is reassigned the value -1, with a warning.
        #' @param data An optional data frame where the variables in
        #'     \code{formula} are to be found. If not found there,
        #'     we search for the variables elsewhere, generally the
        #'     calling environment.
        #' @param likfun An object inheriting from class
        #'     \code{\link{LikelihoodFunction}}, or a generator for such a
        #'     class. This is the likelihood function for the GP model.
        #'     The default is \code{\link{LikLogis}}.
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
                              likfun = LikLogis,      ## likelihood function
                              meanfun = MeanZero,     ## prior mean function
                              covfun = CovSEard,      ## prior cov  function
                              optimize = FALSE,       ## auto optimize hypers?
                              force_all_levs = FALSE, ## keep ALL factor levels?
                              ...) {
            ## Use superclass initializer, making sure to *not* pass "data"
            ## if it is not present.
            if ( missing(data) ) {
                super$initialize(formula, likfun = likfun,
                                 meanfun = meanfun, covfun = covfun,
                                 optimize = optimize,
                                 force_all_levs = force_all_levs)
            } else {
                super$initialize(formula, data, likfun, meanfun, covfun,
                                 optimize, force_all_levs)
            }
            ## Ensure outcomes are in -1, 1
            if ( !all(self$y %in% c(-1, 1)) ) {
                vals <- unique(self$y)
                if ( length(vals) == 2 ) {
                    low_idx  <- which(self$y == min(vals))
                    high_idx <- which(self$y == max(vals))
                    self$y[low_idx]  <- -1
                    self$y[high_idx] <-  1
                    warning("Outcomes automatically translated to {-1, 1}.")
                } else {
                    stop("GPCLA requires dichotomous outcomes.")
                }
            }
            self$train(...)
            if ( optimize ) {
                self$optimize()
            }
        },
        ## Methods
        #' Train the GP model, providing a characterization of the posterior
        #' of the function of interest at the input values given in the
        #' training data.
        #' @details
        #' The implementation here follows very closely Algorithm 3.1 in
        #' Rasmussen and Williams (2006). The main difference is that like
        #' the MATLAB/Octave software GPML, we use successive over-relaxation
        #' in the Newton steps, with an adaptive relaxation factor, though the
        #' implementation differs somewhat from GPML.
        #' @param ... Additional arguments affecting the inference calculations
        #'     (usused).
        train = function(...) {
            ## First a few helper functions:
            ## Psi calculates the negative unnormalized posterior and
            ## allows direct manipulation of alpha in this function itself
            Psi <- function(step_size, alpha, delta, K, m, y) {
                a <- alpha + step_size * delta
                f <- K %*% a + m
                return(c(0.5*crossprod(a,f-m)-sum(plogis(y*f, log.p = TRUE))))
            }
            ## Gets the optimal step size
            step_size <- function(alpha, delta, K, m, y) {
                return(optimize(Psi, c(0, 2), alpha, delta, K, m, y)$minimum)
            }
            ## Now set bookkeeping variables & initialize variables of interest
            n  <- length(self$y)                    ## Number of observations
            m  <- self$meanfun$mean(self$X)         ## Mean function evaluation
            K  <- self$covfun$cov(self$X)           ## Prior covariance
            j  <- rep(1, n)                         ## one vec (for pred probs)
            I  <- diag(n)                           ## An n x n identity matrix
            f  <- m                                 ## Start f at prior mean
            alpha <- rep(0, n)                      ## \nabla log p(y | f)
            o0 <- Inf                               ## Old objective fun eval
            o1 <- Psi(0, alpha, 0, K, m, self$y)    ## Current obj fun eval
            i  <- 0                                 ## Iterations completed
            while ( abs(o1-o0) > 1e-6 & i < 20 ) {  ## Start Newton
                i  <- i  + 1                        ## Update iters completed
                o0 <- o1                            ## Store old obj fun val
                p  <- exp(self$likfun$lp(j, f))     ## Pred probs (MAP)
                W  <- diag(-self$likfun$f_derivative(self$y, f, order = 2))
                sW <- sqrt(W)                       ## Square root of W
                M  <- sW %*% K                      ## We reuse sW %*% K
                L  <- t(chol(I + M %*% sW))         ## LL^T = B = I + sW K sW
                g  <- self$likfun$f_derivative(self$y,f) ## likelihood gradient
                b  <- W %*% (f-m) + g               ## W (f-m) + g
                a  <- b - sW %*% (L %//% (M %*% b)) ## Regular Newtown step
                d  <- c(a - alpha)                  ## Diff btw Newton & alpha
                s  <- step_size(alpha,d,K,m,self$y) ## Find optimal step size
                alpha <- alpha + s * d              ## Update alpha
                o1 <- Psi(0, alpha, 0, K, m,self$y) ## Update obj fun value
                f  <- c(K %*% alpha + m)            ## Update f
            }                                       ## End Newton
            ## Store quantities of future interest
            self$L <- L
            self$alpha <- alpha
            self$sW <- sqrt(diag(-self$likfun$f_derivative(self$y,f,order = 2)))
            self$post_mean <- f
            v <- solve(self$L, self$sW %*% K)
            self$post_cov  <- K - crossprod(v)
            self$prior_mean <- m
            invisible(self)
        },
        #' @description
        #' Characterize the posterior predictive distribution of the function
        #' of interest at new test points.
        #' @details
        #' The implementation here follows very closely Algorithm 3.2 in
        #' Rasmussen and Williams (2006).
        #' @param newdata A data frame containing the data for the new test
        #'     points
        #' @param ... Additional arguments affecting the predictions produced
        predict = function(newdata, ...) {
            X     <- self$X
            Xstar <- construct_Xstar(self, newdata)
            Kss   <- self$covfun$cov(Xstar)
            Ks    <- self$covfun$cov(Xstar, X)
            fmean <- c(self$meanfun$mean(Xstar) + Ks %*% self$alpha)
            fcov  <- solve(self$L, tcrossprod(self$sW, Ks))
            fcov  <- Kss - crossprod(fcov)
            return(
                list(
                    fmean = fmean,
                    fcov  = fcov
                )
            )
        },
        #' @description
        #' Caclulate the negative log marginal likelihood of the GP model.
        #' @param ... Additional arguments affecting the calculation
        nlml = function(...) {
            fhat <- self$post_mean
            kern <- crossprod(fhat - self$prior_mean, self$alpha)
            ll   <- sum(self$likfun$lp(self$y, fhat)) ## log likelihood
            hldB <- sum(log(diag(self$L)))            ## 1/2 log(det(B))
            return(c(0.5 * kern + hldB - ll))
        },
        #' @description
        #' Caclulate the gradient of the negative log marginal likelihood of
        #' the GP model with respect to the hyperparameters of the mean,
        #' covariance, and likelihood functions.
        #' @details
        #' The implementation here follows very closely Algorithm 5.1 in
        #' Rasmussen and Williams (2006).
        #' @param ... Additional arguments affecting the calculation
        dnlml = function(...) {
            ## extract some data members for convenience
            a    <- self$alpha ## \nabla lp, lp=\log p\left(y\mid\hat{f}\right)
            sW   <- self$sW    ## W = -\nabla \nabla lp; sW = W^{1/2}
            L    <- self$L     ## L L^T = B = I + sW K sW
            fhat <- self$post_mean
            ## on to the calculations
            R   <- sW %*% (L %//% sW)      ## = sW (I + sW K sW)^{-1} sW
            K   <- self$covfun$cov(self$X) ## prior covariance
            KR  <- K %*% R                 ## we reuse this a lot
            dl1 <- self$likfun$f_derivative(self$y, fhat, order = 1)
            dl3 <- self$likfun$f_derivative(self$y, fhat, order = 3)
            C   <- solve(L, sW %*% K) ## a step to df = deriv of nlml wrt f, =
            df  <- ((diag(K)-colSums(C^2))/2)*dl3 ## 1/2 diag(K^{-1}W)^{-1}) dl3
            res <- list()
            res[["mean"]] <- numeric(length(self$meanfun$hypers))
            for ( i in seq_along(res[["mean"]]) ) {
                d <- self$meanfun$parameter_derivative(self$X, param = i)
                res[["mean"]][i] <- -crossprod(a, d) - crossprod(df, d - KR%*%d)
            }
            res[["cov"]]  <- numeric(length(self$covfun$hypers))
            for ( i in seq_along(res[["cov"]]) ) {
                d <- self$covfun$parameter_derivative(self$X, param = i, K = K)
                e <- c((sum(colSums(t(R) * d)) - crossprod(a, d) %*% a) / 2)
                b <- d %*% dl1 ## temp calculation
                res[["cov"]][i] <- e - crossprod(df, (b - KR %*% b))
            }
            ## FIXME: At some point I'll need to implement partial derivatives
            ##        of nlml wrt likelihood hypers, but I'm not worrying about
            ##        it right now as all the likelihoods I'd use with binary
            ##        classifications do not have likelihood hypers
            # res[["lik"]]
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
        #' calculated for all predictors).
        #'
        #' A user can request marginal effects on the probability scale,
        #' i.e. change in probability of positive outcome,
        #' by specifying \code{type = "response"},
        #' or marginal effects on the latent function scale,
        #' i.e. change in \eqn{f},
        #' by specifying \code{type = "link"}.
        #' The latter are sometimes called "partial effects."
        #' The partial effects can be calculated directly,
        #' but the marginal effects must be simulated.
        #'
        #' If a predictor is binary, instead the distribution of the difference
        #' in \eqn{f} between the predictor taking the binary variable's
        #' highest value and its lowest value is calculated.
        #' (Or, if \code{type = "response"}, the distribution in the
        #' difference of a positive response is simulated).
        #' If a predictor is categorical, a similar calculation is
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
        #' marginal effects; it calculates the mean, lower bound, and upper
        #' bound (the width of the confidence interval is specifiable) of the
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
        #' @param type A character vector of length one; should be one of
        #'     "link" or "response". If "link", marginal effects are calculated
        #'     on the latent function scale (i.e. change in \eqn{f}). If
        #'     "reponse", marginal effects are calculated on the probability
        #'     scale (i.e. change in probability of positive response).
        #'     The default is "link".
        #' @param M The number of marginal effect draws to take if
        #'     \code{type == "response"}. The default is 1000.
        #' @param ... Other arguments affecting the calculations (unused)
        margins = function(variables = NULL, base_categories = NULL,
                           differences = NULL, indices = NULL, ci = 0.95,
                           force = FALSE, type = "link", M = 1000, ...) {
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
            ## If type == "response", we need the lower Cholesky of K
            if ( type == "response" ) {
                LK <- try(t(chol(K)), silent = TRUE)
                ## If chol() fails, try regularization
                if ( inherits(LK, "try-error") ) {
                    LK <- t(chol(K + diag(1e-3, nrow = n)))
                }
            }
            ## Get effects for each variable
            for ( d in variables ) {
                x   <- self$X[ , d]
                if ( type == "link" ) { ## effects on latent function scale
                    if ( d %in% cnames ) { ## if this is a continuous variable
                        ## Do partial derivatives of prior mean and covariance
                        i   <- which(colnames(self$X) == d)
                        md  <- self$meanfun$input_derivative(self$X,
                                                             dimension = i)
                        Kd  <- self$covfun$input_derivative(self$X,
                                                            dimension = i,
                                                            K = K)
                        Kdd <- self$covfun$input_derivative(self$X,
                                                            dimension = i,
                                                            K = K,
                                                            order = 2)
                        ## Get mean and variance of derivative of f wrt d
                        E   <- c(md + Kd %*% self$alpha)
                        v   <- solve(self$L, self$sW %*% Kd)
                        V   <- Kdd - crossprod(v)
                        rownames(V) <- colnames(V) <- rownames(self$X)
                    } else if ( d %in% bnames ) { ## if it's a binary variable
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
                        v       <- solve(self$L, tcrossprod(self$sW, Ks))
                        Kp      <- Kss - crossprod(v)
                        nn      <- nrow(Xs)
                        i1      <- 1:n
                        i2      <- (n+1):nn
                        V       <- Kp[i1, i1]+Kp[i2, i1]+Kp[i1, i2]+Kp[i2, i2]
                        rownames(V) <- colnames(V) <- rownames(self$X)
                        E       <- c(fs[i1] - fs[i2])
                    } else { ## if this is a categorical variable
                        ## Figure out the baseline category
                        b0 <- names(factors)[which(sapply(names(factors),
                                                          grepl, d))]
                        if ( !is.null(base_categories)
                             & b0 %in% names(base_categories)
                             & isTRUE(base_categories[[b0]] %in% factors[[b0]])
                             ) {
                            b <- paste0(b0, base_categories[[b0]])
                        } else {
                            b <- paste0(b0, factors[[b0]][1])
                        }
                        ## If this *is* the baseline category, skip it
                        if ( d == b ) {
                            next()
                        }
                        ## Get mean and variance of f(d) - f(base_cat)
                        ## First create copy of X where all observations are d
                        X1      <- self$X
                        X1[,paste0(b0, factors[[b0]])] <- 0
                        X1[,d]  <- 1
                        ## Then create copy of X where all observations are b
                        X0      <- self$X
                        X0[,paste0(b0, factors[[b0]])] <- 0
                        X0[,b]  <- 1
                        ## Now get the posterior covariance of those obs
                        Xs      <- rbind(X1, X0)
                        Ks      <- self$covfun$cov(Xs, self$X)
                        Kss     <- self$covfun$cov(Xs)
                        v       <- solve(self$L, tcrossprod(self$sW, Ks))
                        Kp      <- Kss - crossprod(v)
                        ## And the posterior mean
                        fs      <- Ks %*% self$alpha + self$meanfun$mean(Xs)
                        ## And then we can get the mean & cov of the diff
                        nn      <- nrow(Xs)
                        i1      <- 1:n
                        i2      <- (n+1):nn
                        V       <- Kp[i1, i1]+Kp[i2, i1]+Kp[i1, i2]+Kp[i2, i2]
                        rownames(V) <- colnames(V) <- rownames(self$X)
                        E       <- c(fs[i1] - fs[i2])
                    }
                    self$marginal_effects[[d]] <- structure(
                        list(
                            mean = E,
                            covariance = V
                        ),
                        type = type
                    )
                } else { ## effects on probability scale
                    if ( d %in% cnames ) { ## if this is a continuous variable
                        ## Do partial derivatives of prior mean and covariance
                        i   <- which(colnames(self$X) == d)
                        md  <- self$meanfun$input_derivative(self$X,
                                                             dimension = i)
                        Kd  <- self$covfun$input_derivative(self$X,
                                                            dimension = i,
                                                            K = K)
                        Kdd <- self$covfun$input_derivative(self$X,
                                                            dimension = i,
                                                            K = K,
                                                            order = 2)
                        ## Draw f values
                        fhat <- self$post_mean
                        fcov <- self$post_cov
                        f <- rmvn(n = M, mu = fhat, Sigma = fcov)
                        ## Draw f_d values conditional on f values
                        fd <- matrix(NA_real_, nrow = M, ncol = n)
                        v  <- solve(LK, t(Kd))
                        Cd <- Kdd - crossprod(v)
                        fd <- rmvn(n = M, Sigma = Cd)
                        for ( iter in 1:M ) {
                            tmp <- LK %//% (f[iter, ] - self$prior_mean)
                            fd[iter,] <- fd[iter, ] + md + Kd %*% tmp
                        }
                        ## Get pi_{id} values for each draw
                        ## FIXME: This is hard-coded for now,
                        ##        but will need to change
                        draws <- plogis(f) * (1 - plogis(f)) * fd
                    } else if ( d %in% bnames ) { ## if it's a binary variable
                        ## Get draws of sigma(f(1)) - sigma(f(0))
                        ## 1st make copies of X where d is always 1 & always 0
                        if ( d %in% names(differences) ) {
                            vals <- differences[[d]]
                        } else {
                            vals <- sort(unique(x))
                        }
                        X1      <- self$X
                        X1[,d]  <- vals[2]
                        X0      <- self$X
                        X0[,d]  <- vals[1]
                        ## Now take M draws from the f posterior given X1 & X0
                        Ks      <- self$covfun$cov(X1, self$X)
                        fs      <- Ks %*% self$alpha + self$meanfun$mean(X1)
                        Kss     <- self$covfun$cov(X1)
                        v       <- solve(self$L, tcrossprod(self$sW, Ks))
                        Kp      <- Kss - crossprod(v)
                        f1      <- mvtnorm::rmvnorm(M, mean = fs, sigma = Kp)
                        Ks      <- self$covfun$cov(X0, self$X)
                        fs      <- Ks %*% self$alpha + self$meanfun$mean(X0)
                        Kss     <- self$covfun$cov(X0)
                        v       <- solve(self$L, tcrossprod(self$sW, Ks))
                        Kp      <- Kss - crossprod(v)
                        f0      <- mvtnorm::rmvnorm(M, mean = fs, sigma = Kp)
                        ## Now push through the sigmoid & store the difference
                        draws   <- plogis(f1) - plogis(f0)
                    } else { ## if this is a categorical variable
                        ## Figure out the baseline category
                        b0 <- names(factors)[which(sapply(names(factors),
                                                          grepl, d))]
                        if ( !is.null(base_categories)
                             & b0 %in% names(base_categories)
                             & isTRUE(base_categories[[b0]] %in% factors[[b0]])
                        ) {
                            b <- paste0(b0, base_categories[[b0]])
                        } else {
                            b <- paste0(b0, factors[[b0]][1])
                        }
                        ## If this *is* the baseline category, skip it
                        if ( d == b ) {
                            next()
                        }
                        ## Get mean and variance of f(d) - f(base_cat)
                        ## First create copy of X where all observations are d
                        X1      <- self$X
                        X1[,paste0(b0, factors[[b0]])] <- 0
                        X1[,d]  <- 1
                        ## Then create copy of X where all observations are b
                        X0      <- self$X
                        X0[,paste0(b0, factors[[b0]])] <- 0
                        X0[,b]  <- 1
                        ## Now take M draws from the f posterior given X1 & X0
                        Ks      <- self$covfun$cov(X1, self$X)
                        fs      <- Ks %*% self$alpha + self$meanfun$mean(X1)
                        Kss     <- self$covfun$cov(X1)
                        v       <- solve(self$L, tcrossprod(self$sW, Ks))
                        Kp      <- Kss - crossprod(v)
                        f1      <- mvtnorm::rmvnorm(M, mean = fs, sigma = Kp)
                        Ks      <- self$covfun$cov(X0, self$X)
                        fs      <- Ks %*% self$alpha + self$meanfun$mean(X0)
                        Kss     <- self$covfun$cov(X0)
                        v       <- solve(self$L, tcrossprod(self$sW, Ks))
                        Kp      <- Kss - crossprod(v)
                        f0      <- mvtnorm::rmvnorm(M, mean = fs, sigma = Kp)
                        ## Now push through the sigmoid & store the difference
                        draws   <- plogis(f1) - plogis(f0)
                    }
                    self$marginal_effects[[d]] <- structure(
                        list(
                            draws = draws
                        ),
                        type = type
                    )
                } ## End variable type conditional
            } ## End for loop over variables
            ## Average the marginal effects over the requested indices
            if ( is.null(indices) ) {
                indices <- 1:length(self$y)
            }
            ames <- sapply(self$marginal_effects, function(x) {
                if ( attr(x, "type") == "link" ) {
                    return( mean(x[["mean"]][indices]) )
                } else {
                    return( mean(colMeans(x[["draws"]][ , indices])) )
                }
            })
            names(ames) <- names(self$marginal_effects)
            res <- data.frame(
                Variable = names(ames),
                Mean = ames,
                LB = 0, UB = 0
            )
            k <- 1/(length(indices)^2)
            ci_low <- (1-ci)/2
            ci_high <- 1-((1-ci)/2)
            for ( i in seq_len(nrow(res)) ) {
                d <- res$Variable[i]
                x <- self$marginal_effects[[d]]
                if ( attr(x, "type") == "link" ) {
                    this_sd <- sqrt(k*sum(x[["covariance"]][indices, indices]))
                    res$LB[i] <- qnorm(ci_low,  mean = ames[d], sd = this_sd)
                    res$UB[i] <- qnorm(ci_high, mean = ames[d], sd = this_sd)
                } else {
                    draws <- rowMeans(x[["draws"]][ , indices])
                    res$LB[i] <- quantile(draws, probs = ci_low)
                    res$UB[i] <- quantile(draws, probs = ci_high)
                }
            }
            self$average_marginal_effects <- res
            invisible(self)
        }
    )
)
