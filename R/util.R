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
