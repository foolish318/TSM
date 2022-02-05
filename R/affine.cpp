// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

/* 
These functions calculate affine loadings for a Gaussian Dynamic Term Structure Model. 
 */

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
List affineLoadingsCpp(vec mu, mat phi, mat Omega, double delta0, vec delta1, vec mats, double dt) {
    int J = mats.n_elem;
    int N = delta1.n_elem;
    mat A = zeros(1, J);
    mat B = zeros(N, J);
    double Atmp = 0;
    vec Btmp = zeros(N);
    mat muprime = mu.t();
    mat phiprime = phi.t();
    int i = 0; // index of mats
    double ip1, ip2;
    for (int n = 1; n <= mats.max(); n++) {
        ip1 = as_scalar(muprime * Btmp);
        ip2 = as_scalar(0.5 * Btmp.t() * Omega * Btmp);
        Atmp = Atmp + ip1 + ip2 - delta0;
        Btmp = phiprime * Btmp - delta1;
        if (n == mats(i)) {
            A(i) = - Atmp/n;
            B(span::all, i) = - Btmp/n;
            i++;
        }
    }
    return List::create(Named("A") = A/dt,
                        Named("B") = B/dt);
}

// [[Rcpp::export]]
List affineLoadingsBonlyCpp(mat phi, vec delta1, vec mats, double dt) {
    // if no intercept (A matrix) is required, this saves some time
    int J = mats.n_elem;
    int N = delta1.n_elem;
    mat B = zeros(N, J);
    vec Btmp = zeros(N);
    mat phiprime = phi.t();
    int i = 0; // index of mats
    for (int n = 1; n <= mats.max(); n++) {
        Btmp = phiprime * Btmp - delta1;
        if (n == mats(i)) {
            B(span::all, i) = - Btmp/n;
            i++;
        }
    }
    return List::create(Named("B") = B/dt);
}
