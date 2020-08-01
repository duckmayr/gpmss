#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat cholsove(const arma::mat& L, const arma::mat& B) {
    using arma::solve;
    using arma::trimatl;
    using arma::trimatu;
    return solve(trimatu(L.t()), solve(trimatl(L), B));
}
