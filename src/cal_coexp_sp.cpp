#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' cal_coexp_sp
//' This function calculates the coexpression patterns between genes
//' in sparse matrix and returns the coexpression matrix.
//' @author Qi Gao
//' @param X Input binarized cell (row) by gene (column) sparse matrix
//' @return Coexpression matrix
//' @export
// [[Rcpp::export]]
arma::mat cal_coexp_sp(arma::sp_mat X, arma::sp_mat X_subcohort){
  int p = X.n_cols;
  int n = X.n_rows;
  arma::vec q(p);
  for(int i = 0; i < p; i++){
    q(i) = mean(X.col(i));
  }
  arma::vec mq = 1 - q;
  arma::mat c = X_subcohort.t() * X_subcohort - q * q.t() * n;
  arma::mat d = sqrt(X_subcohort.n_rows * q * q.t() % (mq * mq.t()));
  
  return(c / d);
}

