context("Covariance Functions")

##### Setup -----
## Setup data; simple 3x2 and 2x2 matrices that are easy to know right answers
X <- matrix(1:6, nrow = 3, ncol = 2)
Z <- matrix(7:0, nrow = 2, ncol = 2)
n <- nrow(X)


##### Isotropic Squared Exponential -----
## Create the CovarianceFunction object
h <- c(1, -1)
k <- CovSEiso$new(hypers = h)
## Test that it's creating covariance matrices correctly
test_that("CovSEiso's covariance method works properly", {
    ## Test cov(X, X)
    K <- k$cov(X)
    ## Diagonal should be scale factor
    expect_equal(diag(K), rep(exp(2 * h[1]), n))
    ## It should be symmetric
    expect_true(isSymmetric(K))
    ## Do a manual check on one of its values
    K_n_1_manual <- exp(2*h[1]) * exp(-0.5 * sum((X[n,]-X[1,])^2) / exp(2*h[2]))
    expect_equal(K[3,1], K_n_1_manual)
    ## Test cov(X, Z)
    Ks <- k$cov(X, Z)
    Ks_1_1_manual <- exp(2*h[1]) * exp(-0.5*sum((X[1,]-Z[1,])^2) / exp(2*h[2]))
    ## Do a manual check on one of its values
    expect_equal(Ks[1,1], Ks_1_1_manual)
    ## Check its dimensions
    expect_equal(dim(Ks), c(nrow(X), nrow(Z)))
    ## Check the check on hypers
    expect_error(k$cov(X, X, 1))
})
## Test the partial derivatives with respect to hypers
test_that("CovSEiso's derivatives with respect to parameters work properly", {
    ## Check partial derivative wrt log(sigma_f)
    d_sf <- k$parameter_derivative(X, X, h, 1)
    K    <- k$cov(X)
    expect_equal(d_sf, 2 * K)
    ## Check partial derivative wrt log(ell)
    d_ell <- k$parameter_derivative(X, X, h, 2)
    expect_equal(diag(d_ell), rep(0, n))
    expect_equal(d_ell[n, 1], K[n, 1] * sum(((X[n, ] - X[1, ])^2)/exp(2*h[2])))
    ## Check the check on hypers
    expect_error(k$parameter_derivative(X, X, 1))
})
## Test the partial derivatives with respect to inputs
test_that("CovSEiso's derivatives with respect to inputs work properly", {
    ## Check partial derivative wrt dim 1 of X
    d_d1 <- k$input_derivative(X, X, h, 1, 1)
    K    <- k$cov(X)
    expect_equal(d_d1[n, 1], K[n, 1] * (X[1, 1] - X[n, 1]) / exp(2*h[2]))
    ## Check cross partial derivative wrt dim 1 of X
    d_d2 <- k$input_derivative(X, X, h, 1, 2)
    ell2 <- exp(2 * h[2])
    expect_true(isSymmetric(d_d2))
    expect_equal(diag(d_d2), rep(exp(2 * h[1]) / ell2, n))
    d <- X[n, 1] - X[1, 1]
    expect_equal(d_d2[n, 1], K[n, 1] * (1/ell2 + (d*-d)/(ell2*ell2)))
    ## Check the check on hypers
    expect_error(k$input_derivative(X, X, 1))
    ## Check the check on derivative order
    expect_error(k$input_derivative(X, X, h, 1, 3))
})


##### Squared Exponential with Automatic Relevance Determination -----
## Create the CovarianceFunction object
h <- c(1, -1, 1)
k <- CovSEard$new(hypers = h)
## Test that it's creating covariance matrices correctly
test_that("CovSEard's covariance method works properly", {
    ## Test cov(X, X)
    K <- k$cov(X)
    ## Diagonal should be scale factor
    expect_equal(diag(K), rep(exp(2 * h[1]), n))
    ## It should be symmetric
    expect_true(isSymmetric(K))
    ## Do a manual check on one of its values
    K_n_1_manual <- exp(2*h[1]) * exp(-0.5*sum((X[n,]-X[1,])^2 / exp(2*h[-1])))
    expect_equal(K[3,1], K_n_1_manual)
    ## Test cov(X, Z)
    Ks <- k$cov(X, Z)
    ## Do a manual check on one of its values
    Ks_1_1_manual <- exp(2*h[1]) * exp(-0.5*sum((X[1,]-Z[1,])^2 / exp(2*h[-1])))
    expect_equal(Ks[1,1], Ks_1_1_manual)
    ## Check its dimensions
    expect_equal(dim(Ks), c(nrow(X), nrow(Z)))
    ## Check the check on hypers
    expect_error(k$cov(X, X, 1))
})
## Test the partial derivatives with respect to hypers
test_that("CovSEard's derivatives with respect to parameters work properly", {
    ## Check partial derivative wrt log(sigma_f)
    d_sf <- k$parameter_derivative(X, X, h, 1)
    K    <- k$cov(X)
    expect_equal(d_sf, 2 * K)
    ## Check partial derivative wrt log(ell)
    d_ell <- k$parameter_derivative(X, X, h, 2)
    expect_equal(diag(d_ell), rep(0, n))
    expect_equal(d_ell[n, 1], K[n, 1] * sum(((X[n, 2]-X[1, 2])^2)/exp(2*h[2])))
    ## Check the check on hypers
    expect_error(k$parameter_derivative(X, X, 1))
})
## Test the partial derivatives with respect to inputs
test_that("CovSEard's derivatives with respect to inputs work properly", {
    ## Check partial derivative wrt dim 1 of X
    d_d1 <- k$input_derivative(X, X, h, 1, 1)
    K    <- k$cov(X)
    expect_equal(d_d1[n, 1], K[n, 1] * (X[1, 1] - X[n, 1]) / exp(2*h[2]))
    ## Check cross partial derivative wrt dim 1 of X
    d_d2 <- k$input_derivative(X, X, h, 1, 2)
    ell2 <- exp(2 * h[2])
    expect_true(isSymmetric(d_d2))
    expect_equal(diag(d_d2), rep(exp(2 * h[1]) / ell2, n))
    d <- X[n, 1] - X[1, 1]
    expect_equal(d_d2[n, 1], K[n, 1] * (1/ell2 + (d*-d)/(ell2*ell2)))
    ## Check the check on hypers
    expect_error(k$input_derivative(X, X, 1))
    ## Check the check on derivative order
    expect_error(k$input_derivative(X, X, h, 1, 3))
})
