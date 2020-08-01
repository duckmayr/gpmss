#include <RcppArmadillo.h>

// [[Rcpp::export(.covSEiso)]]
arma::mat covSEiso(const arma::mat& X, const arma::mat& Z, const arma::vec& h) {
    arma::uword n = X.n_rows;
    arma::uword m = Z.n_rows;
    arma::uword p = X.n_cols;
    double sf2    = std::exp(2.0 * h[0]);
    double l      = -0.5 / std::exp(2.0 * h[1]);
    arma::mat res(n, m);
    for ( arma::uword j = 0; j < m; ++j ) {
        for ( arma::uword i = 0; i < n; ++i ) {
            double elem = 0.0;
            for ( arma::uword k = 0; k < p; ++k ) {
                double diff = X(i, k) - Z(j, k);
                elem += (diff * diff);
            }
            res(i, j) = sf2 * std::exp(l * elem);
        }
    }
    return res;
}
