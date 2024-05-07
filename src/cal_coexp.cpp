#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' cal_coexp
//' This function calculates the coexpression patterns between genes
//' and returns the coexpression matrix.
//' @author Qi Gao
//' @param X Input binarized cell (row) by gene (column) matrix
//' @return Coexpression matrix
//' @export
// [[Rcpp::export]]
arma::mat cal_coexp(arma::mat X, arma::mat X_subcohort){
  int p = X_subcohort.n_cols;
  int n = X_subcohort.n_rows;
  arma::vec q(p);
  for(int i = 0; i < p; i++){
    q(i) = mean(X_subcohort.col(i));
  }
  arma::vec mq = 1 - q;
  arma::mat c = X_subcohort.t() * X_subcohort - q * q.t() * X_subcohort.n_rows;
  arma::mat d = sqrt(X_subcohort.n_rows * q * q.t() % (mq * mq.t()));
  
  return(c / d);
}
