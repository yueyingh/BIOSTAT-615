#include "utils.h"
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List admm_genlasso(const arma::mat &A, const arma::colvec &b,
                         const arma::mat &D, const arma::mat &M,
                         const double lambda,
                         const double reltol, const double abstol,
                         const int maxiter, const double rho) {
  // 1. get parameters
  //   const int m = A.n_rows;
  const int n = A.n_cols;

  // 2. set ready
  arma::colvec x(n, fill::randn);
  x /= 10.0;
  arma::colvec z(D * x);
  arma::colvec u_prev(D * x - z);
  arma::colvec u(D * x - z);
  arma::colvec q(n, fill::zeros);
  arma::colvec zold(z);
  arma::colvec x_hat(n, fill::zeros);

  // 3. precompute static variables for x-update and factorization
  arma::mat Atb = A.t() * b;
  arma::mat U = genlasso_factor(A, rho, M); // returns upper
  arma::mat L = U.t();

  // 4. iteration
  arma::vec h_objval(maxiter, fill::zeros);
  arma::vec h_r_norm(maxiter, fill::zeros);
  arma::vec h_s_norm(maxiter, fill::zeros);
  arma::vec h_eps_pri(maxiter, fill::zeros);
  arma::vec h_eps_dual(maxiter, fill::zeros);

  double sqrtn = std::sqrt(static_cast<float>(n));
  int k;
  for (k = 0; k < maxiter; k++) {
    // 4-1. update 'x'
    q = Atb-D.t()*(2*u-u_prev)+rho/2*M*x+rho/2*M.t()*x; // temporary value
    x = solve(trimatu(U), solve(trimatl(L), q));
    //        if (m >= n){
    //            x = solve(trimatu(U),solve(trimatl(L),q));
    //        } else {
    //            x = q/rho -
    //            (A.t()*solve(trimatu(U),solve(trimatl(L),A*q)))/rho2;
    //        }

    // 4-2. update 'z'
    zold = z;
    z = genlasso_shrinkage(D * x + u/rho, lambda / rho);

    // 4-3. update 'u'
    u_prev = u
    u = u + rho * (D * x - z);

    // 4-3. dianostics, reporting
    h_objval(k) = genlasso_objective(A, b, D, lambda, x, z);
    h_r_norm(k) = arma::norm(D * x - z);
    h_s_norm(k) = arma::norm(-rho * (z - zold));
    if (norm(x) > norm(-z)) {
      h_eps_pri(k) = sqrtn * abstol + reltol * norm(x);
    } else {
      h_eps_pri(k) = sqrtn * abstol + reltol * norm(-z);
    }
    h_eps_dual(k) = sqrtn * abstol + reltol * norm(rho * u);

    // 4-4. termination
    if ((h_r_norm(k) < h_eps_pri(k)) && (h_s_norm(k) < h_eps_dual(k))) {
      break;
    }
  }

  // 5. report results
  List output;
  output["x"] = x;             // coefficient function
  output["objval"] = h_objval; // |x|_1
  output["k"] = k;             // number of iterations
  output["r_norm"] = h_r_norm;
  output["s_norm"] = h_s_norm;
  output["eps_pri"] = h_eps_pri;
  output["eps_dual"] = h_eps_dual;
  return (output);
}