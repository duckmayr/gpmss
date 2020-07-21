context("Likelihood Functions")

##### Setup -----
## Setup data
y <- c(1, -1,  1)
f <- c(5, -2, -1)

##### Normal Likelihood -----
lik <- LikGauss$new(1)
test_that("LikGauss's lp() method works properly", {
    lp <- dnorm(y, mean = f, sd = exp(1), log = TRUE)
    expect_equal(lik$lp(y, f), lp)
    expect_error(lik$lp(y, 1))
})
test_that("LikGauss's f_derivative() method works properly", {
    expect_null(lik$f_derivative(y, f))
})

##### Logistic Likelihood -----
lik <- LikLogis$new()
test_that("LikLogis's lp() method works properly", {
    lp <- plogis(y * f, log = TRUE)
    expect_equal(lik$lp(y, f), lp)
    expect_error(lik$lp(y, 1))
})
test_that("LikLogis's f_derivative() method works properly", {
    p    <- plogis(f)
    d1lp <- ( (y+1) / 2 ) - p
    d2lp <- -p * (1 - p)
    d3lp <- 2 * d2lp * (0.5 - p)
    expect_equal(lik$f_derivative(y, f, order = 1), d1lp)
    expect_equal(lik$f_derivative(y, f, order = 2), d2lp)
    expect_equal(lik$f_derivative(y, f, order = 3), d3lp)
    expect_error(lik$f_derivative(y, f, order = 4))
    expect_error(lik$f_derivative(y, 1))
})
