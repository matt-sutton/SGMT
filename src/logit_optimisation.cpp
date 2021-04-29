// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>
using namespace Rcpp;
using namespace Numer;
using namespace Eigen;


typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

class LogisticReg: public MFuncGrad
{
private:
    const MapMat X;
    const MapVec Y;
    const MapVec UmZ;
    const double Rho;
    const int n;
    Eigen::VectorXd xbeta;  // contains X*beta
    Eigen::VectorXd prob;   // contains log(1+exp(X*beta)) and 1/(1+exp(-X*beta))
public:
    LogisticReg(const MapMat x_, const MapVec y_, const MapVec umz_, const double rho_) :
    X(x_), Y(y_), UmZ(umz_), Rho(rho_), n(X.rows()), xbeta(n), prob(n) {}

    double f_grad(Constvec& beta, Refvec grad)
    {
        // Negative log likelihood
        //   sum(log(1 + exp(X * beta))) - y' * X * beta

        xbeta.noalias() = X * beta;
        const double yxbeta = Y.dot(xbeta);
        // Calculate log(1 + exp(X * beta)), avoiding overflow
        for(int i = 0; i < n; i++)
            prob[i] = R::log1pexp(xbeta[i]); // should find a better solution for this

        const double f = prob.sum() - yxbeta + (Rho/2)*( beta + UmZ ).squaredNorm();

        // Gradient
        //   X' * (p - y), p = exp(X * beta) / (1 + exp(X * beta))

        // exp(X * beta) => p
        prob = (xbeta - prob).array().exp();
        grad.noalias() = X.transpose() * (prob - Y) + Rho*( beta + UmZ );

        return f;
    }
};

// [[Rcpp::export]]
Rcpp::List b_update(Rcpp::NumericMatrix x, Rcpp::NumericVector y, Rcpp::NumericVector umz, double rho, Rcpp::NumericVector beta_init)
{
    const MapMat xx = Rcpp::as<MapMat>(x);
    const MapVec yy = Rcpp::as<MapVec>(y);
    const MapVec uumzz = Rcpp::as<MapVec>(umz);
    const double rrho = rho;
    // Negative log likelihood
    LogisticReg nll(xx, yy, uumzz, rrho);
    // Initial guess
    Rcpp::NumericVector b = Rcpp::clone(beta_init);
    MapVec beta(b.begin(), b.length());

    double fopt;
    int status = optim_lbfgs(nll, beta, fopt, 100, 1e-8, 1e-5);
    if(status < 0)
        Rcpp::stop("fail to converge");

    return Rcpp::List::create(
        Rcpp::Named("coefficients") = b,
        Rcpp::Named("converged")    = (status >= 0));
}
