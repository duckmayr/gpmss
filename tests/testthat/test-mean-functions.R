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
mu1   <- MeanLinear$new()
mu2   <- MeanLinear$new(intercept = FALSE)
alpha <- 1
beta  <- rep(c(1, -1), length.out = ncol(X))
hyp   <- c(alpha, beta)
test_that("MeanLinear's mean() method works properly", {
    expect_equal(mu1$mean(X, hyp), c(alpha + X %*% beta))
})
test_that("MeanLinear's parameter_derivative() method works properly", {
    expect_equal(mu1$parameter_derivative(X, hyp, 1), rep(0, nrow(X)))
    expect_equal(mu1$parameter_derivative(X, hyp, 2), X[, 1])
    expect_equal(mu2$parameter_derivative(X, beta, 1), X[, 1])
    expect_error(mu1$parameter_derivative(X, hyp, ncol(X)+2))
})
test_that("MeanLinear's input_derivative() method works properly", {
    expect_equal(mu1$input_derivative(X, hyp, 1), rep(beta[1], n))
    expect_error(mu1$input_derivative(X, hyp, ncol(X)+1))
})
