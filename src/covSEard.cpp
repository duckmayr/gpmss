#include <RcppArmadillo.h>

// [[Rcpp::export(.covSEard)]]
arma::mat covSEard(const arma::mat& X, const arma::mat& Z, const arma::vec& h) {
    arma::uword n = X.n_rows;
    arma::uword m = Z.n_rows;
    arma::uword p = X.n_cols;
    double sf2    = std::exp(2.0 * h[0]);
    arma::vec l   = 1 / arma::exp(2.0 * h.tail(h.n_elem - 1));
    arma::mat res(n, m);
    for ( arma::uword j = 0; j < m; ++j ) {
        for ( arma::uword i = 0; i < n; ++i ) {
            double elem = 0.0;
            for ( arma::uword k = 0; k < p; ++k ) {
                double diff = X(i, k) - Z(j, k);
                elem += (l[k] * diff * diff);
            }
            res(i, j) = sf2 * std::exp(-0.5 * elem);
        }
    }
    return res;
}
