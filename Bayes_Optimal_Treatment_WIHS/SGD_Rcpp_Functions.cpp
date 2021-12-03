//[[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>


using namespace std;
using namespace arma;
using namespace Rcpp; 


// [[Rcpp::export]]
mat rmvn_rcpp(const int n, vec& mean, mat& sigma){
  
  // randomly generate samples from multivariate normal distribution
  // args: n: number of data 
  //      mean: row vector mean, sigma: covaraicne matrix  
  // returns: out: random samples from multivariate normal distribution
  
  int k = sigma.n_cols; // dimension of the multivariate normal distribution
  mat z = randn(n, k);
  mat out = repmat(mean,1,n).t()+z*chol(sigma);
  return(out);
}



// [[Rcpp::export]]
mat auto_regressive_rcpp(const int n, const double sigma2, const double rho, vec& t){
  
  // construct auto regressive covariance matrix 
  // args: n: dimension of the matrix
  //       sigma2: variance, rho: correlation, t: time 
  // returns: Cov: auto regressive covariance matrix 
  
  mat Cov(n, n);
  Cov.fill(sigma2);
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      Cov(i,j) *= pow(rho, abs(t(i)-t(j)));
    }
  }
  return(Cov);
} 



// [[Rcpp::export]]
double noncentral_hypergeo_integrand_rcpp(const double x, const int Nk_k, vec& a2k, vec& tauk, 
        vec& omegak, vec& omegak_dev, const double omegak_sum, const double omegak_devsum, const double ym){
  
  // integrand function for calculating (derivatives) noncentral hypergeometric distribution 
  // args: x: integrand is a function of x 
  //       Nk_k: number of single ART drugs in k-th drug class
  //       a2k: second level tree representation for drug combination in k-th drug class
  //       tauk, omegak, omegak_dev, omega_sum, omgeak_devsum: parameters in noncentral hypergeometric distribution
  //       ym: m-th state value y 
  // returns: res: evaluation of integrand function at x
  
  double res = 0;
  double tauk_dev=0, x_tauk_prod = 1, eps = 1e-10;
  vec x_tauk(Nk_k); 
  
  for (int nk=0; nk<Nk_k; nk++){
    x_tauk[nk] = 1 - pow(x, tauk[nk]);
    x_tauk_prod *= pow(x_tauk[nk], a2k[nk]);
  }
  
  for (int nk_tilde=0; nk_tilde<Nk_k; nk_tilde++){
    // derivative of tauk[nk_tilde] on theta2k[nk]
    tauk_dev = (omegak_dev[nk_tilde]*omegak_sum - omegak[nk_tilde]*omegak_devsum)/pow(omegak_sum,2)*ym;
    // sum the integrand over nk_tilde from 1 to Nk[k] 
    if (x_tauk[nk_tilde] > eps){
      res += tauk_dev * a2k[nk_tilde]*pow(x_tauk[nk_tilde], a2k[nk_tilde]-1) *
        -pow(x, tauk[nk_tilde])*log(x) * x_tauk_prod/pow(x_tauk[nk_tilde], a2k[nk_tilde]);
    }
  }
  
  return(res);
}



// [[Rcpp::export]]
double noncentral_hypergeo_integrate_rcpp(const int num, const double lower, const double upper, const int Nk_k, vec& a2k, 
          vec& tauk, vec& omegak, vec& omegak_dev, const double omegak_sum, const double omegak_devsum, const double ym){
  
  // integrate function for calculating (derivatives) noncentral hypergeometric distribution 
  // args: num: number of bins for integrate approximation
  //       lower, upper: lower and upper bounds for integration
  //       Nk_k: number of single ART drugs in k-th drug class
  //       a2k: second level tree representation for drug combination in k-th drug class
  //       tauk, omegak, omegak_dev, omega_sum, omgeak_devsum: parameters in noncentral hypergeometric distribution
  //       ym: m-th state value y 
  // returns: res: integration of (derivatives) noncentral hypergeometric distribution density over [lower, upper]
  
  double res = 0; 
  double bin = (upper - lower)/num;
  double a = lower, b = lower + bin;
  
  for (int i=0; i<num; i++){
    res += (b-a) * noncentral_hypergeo_integrand_rcpp((b+a)/2, Nk_k, a2k, tauk, omegak, omegak_dev, omegak_sum, omegak_devsum, ym);
    a += bin; b += bin;
  }
  
  return(res);
}



// [[Rcpp::export]]
double drug_toxic_integrand_rcpp(const double x, const int nk, const int Ji, 
                                 mat& Dh, vec& ti, const double upper){
  
  // integrand function for calculating drug_toxic
  // args: x: integrand is a function of x 
  //       nk: the nk-th single ART drug
  //       Ji: number of observations for a single subject
  //       Dh: drug history binary indicator matrix 
  //       ti: (transformed) time vector for a single subject
  //       upper: upper bound for integration
  // returns: res: evaluation of integrand function at x 
  
  double res = 0;
  double drug_usage = 0;
  
  for (int j=0; j<Ji; j++){
    if (j == 0){
      if (x < ti[j]){
        drug_usage = Dh(j, nk);
      }
    }else{
      if (x > ti[j-1] && x < ti[j]){
        drug_usage = Dh(j, nk);
      }
    }
  }
  
  res = drug_usage * exp(-(upper - x));
  return(res);
}



// [[Rcpp::export]]
vec drug_toxic_integrate_rcpp(const int num, const double lower, const double upper, 
                              const int Nk_sum, const int Ji, mat& Dh, vec& ti){
  
  // integrate function for calculating drug_toxic
  // args: num: number of bins for integrate approximation
  //       lower, upper: lower and upper bounds for integration
  //       Nk_sum: number of single ART drugs
  //       Ji: number of observations for a single subject
  //       Dh: drug history binary indicator matrix 
  //       ti: (transformed) time vector for a single subject
  // returns: drug_toxic: vector of drug toxicity for each single ART drug
  
  vec drug_toxic(Nk_sum); drug_toxic.fill(0);
  
  for (int nk=0; nk<Nk_sum; nk++){
    
    // calculate drug toxicity for each single ART drug
    double res = 0; 
    double bin = (upper - lower)/num;
    double a = lower, b = lower + bin;
    
    for (int i=0; i<num; i++){
      res += (b-a) * drug_toxic_integrand_rcpp((b+a)/2, nk, Ji, Dh, ti, upper);
      a += bin; b += bin;
    }
    drug_toxic[nk] = res;
  }
  
  return(drug_toxic);
}