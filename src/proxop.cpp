#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat Lasso(arma::mat x, double lambda) {
  return arma::sign(x)%arma::max(arma::abs(x) - lambda, arma::zeros(arma::size(x)));
}

// [[Rcpp::export]]
arma::mat L21(arma::mat x, double lambda) {
  int p = x.n_rows;
  int q = x.n_cols;
  arma::mat out(p, q);
  for (int g = 0; g < p; g++){
    double nrm = std::max(arma::norm(x.row(g),"fro"), 1e-10);
    out.row(g) = x.row(g)*std::max(1 - lambda/nrm, 0.0);
  }
  return out;
}
// [[Rcpp::export]]
arma::mat L21_w(arma::mat x, arma::vec lambda) {
  int p = x.n_rows;
  int q = x.n_cols;
  arma::mat out(p, q);
  for (int g = 0; g < p; g++){
    double nrm = std::max(arma::norm(x.row(g),"fro"), 1e-10);
    out.row(g) = x.row(g)*std::max(1 - lambda(g)/nrm, 0.0);
  }
  return out;
}
// [[Rcpp::export]]
arma::mat gLasso(arma::mat x, double lambda) {
  double nrm = std::max(arma::norm(x,"fro"), 1e-10);
  return x*std::max(1 - lambda/nrm, 0.0);
}

// [[Rcpp::export]]
arma::mat proxop(arma::mat x, arma::vec groups, double lambda, double alpha) {
  arma::vec gunique = unique(groups);
  double ngroup = gunique.size();
  int p = x.n_rows;
  int q = x.n_cols;
  arma::mat out(p, q);
  double glambda = lambda*(1-alpha);
  double llambda = lambda*alpha;

  for (int g = 0; g < ngroup; g++){
    arma::uvec gind = arma::find(gunique(g) == groups);
    out.rows(gind) = gLasso(L21(x.rows(gind),llambda), glambda*std::sqrt(gind.size()));
  }
  return out;
}

// [[Rcpp::export]]
arma::mat weighted_proxop(arma::mat x, arma::vec groups, double lambda, double alpha, arma::vec grp_weights, arma::vec ind_weights) {
  arma::vec gunique = unique(groups);
  double ngroup = gunique.size();
  int p = x.n_rows;
  int q = x.n_cols;
  arma::mat out(p, q);
  double glambda = lambda*(1-alpha);
  double llambda = lambda*alpha;

  for (int g = 0; g < ngroup; g++){
    arma::uvec gind = arma::find(gunique(g) == groups);
    out.rows(gind) = gLasso(L21_w(x.rows(gind),ind_weights.elem(gind)*llambda), glambda*grp_weights(g));
  }
  return out;
}
