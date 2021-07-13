## Use double left division to solve A^{-1} B = X for X, where L L^T = A
## (unfortunately can't use \\ due to escape character issues)
`%//%` <- function(L, B) {
    return(solve(t(L), solve(L, B)))
}

## The following helper function takes a formula and appends -1
## It does that by converting the formula to a character vector of length one
## (named char_formula), pasting "- 1" on the end, then converting that
## character to a formula.
## It will be useful for getting one-hot encoding for factors
## without the user having to worry about omitting the intercept.
convert_formula <- function(formula_in) {
    char_formula <- paste(trimws(format(formula_in)), collapse = " ")
    return(stats::formula(paste(char_formula, "- 1")))
}

## This helper function generates a design matrix for test data given a
## previously generated model (useful for prediction)
#' @importFrom stats terms
construct_Xstar <- function(model, newdata) {
    Xstar <- data.frame(row.names = seq_len(nrow(newdata)))
    for ( i in attr(terms(model), "term.labels") ) {
        if ( i %in% names(attr(model$X, "contrasts")) ) {
            cntr <- attr(model$X,"contrasts")[[i]]
            Xstar[ , paste0(i, rownames(cntr))] <- cntr[newdata[[i]], ]
        } else {
            Xstar[ , i] <- newdata[[i]]
        }
    }
    return(as.matrix(Xstar))
}

## This function gets the (upper) Cholesky decomposition of a matrix that's
## PSD by construction, but which may have very small eigenvalues;
## we try to decomp and if it fails, we try again after regularizing
safe_chol <- function(M) {
    R <- try(chol(M), silent = TRUE)
    r <- 1e-6
    while ( inherits(R, "try-error") & r < 1 ) {
        R <- try(chol(M + diag(r, nrow = nrow(M))), silent = TRUE)
        r <- r * 10
    }
    if ( inherits(R, "try-error") ) {
        stop("Cholesky decomp failed; matrix was not numerically PSD.")
    }
    return(R)
}

## This is a custom multivariate normal generating function;
## I set this up merely to use the safe_chol() function above
rmvn <- function(n, mu = rep(0, nrow(Sigma)), Sigma) {
    m <- ncol(Sigma)
    R <- safe_chol(Sigma)
    U <- matrix(stats::rnorm(n*m), nrow = n)
    res <- U %*% R
    if ( !all(mu == 0) ) {
        for ( i in 1:n ) {
            res[i, ] <- res[i, ] + mu
        }
    }
    return(res)
}
