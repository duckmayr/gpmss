context("Mean Functions")

##### Setup -----
## Setup data; a simple 3x2 matrix for which it's easy to know right answers
X <- matrix(1:6, nrow = 3, ncol = 2)
n <- nrow(X)

##### Zero Mean -----
mu <- MeanZero$new()
test_that("MeanZero's mean() method works properly", {
    expect_equal(mu$mean(X), rep(0, n))
})
test_that("MeanZero's parameter_derivative() method works properly", {
    expect_equal(mu$parameter_derivative(X), numeric())
})
test_that("MeanZero's input_derivative() method works properly", {
    expect_equal(mu$input_derivative(X), rep(0, n))
})

##### Linear Mean -----
mu   <- MeanLinear$new()
beta <- rep(c(1, -1), length.out = ncol(X))
test_that("MeanLinear's mean() method works properly", {
    expect_equal(mu$mean(X, beta), X %*% beta)
})
test_that("MeanLinear's parameter_derivative() method works properly", {
    expect_equal(mu$parameter_derivative(X, beta, 1), X[, 1])
    expect_error(mu$parameter_derivative(X, beta, ncol(X)+1))
})
test_that("MeanLinear's input_derivative() method works properly", {
    expect_equal(mu$input_derivative(X, beta, 1), rep(beta[1], n))
    expect_error(mu$input_derivative(X, beta, ncol(X)+1))
})
