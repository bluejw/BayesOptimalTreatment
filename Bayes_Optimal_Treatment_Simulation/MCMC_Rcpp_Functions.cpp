//[[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp; 


const double log2pi = log(2.0*M_PI);


// [[Rcpp::export]]
double dmvn_rcpp(rowvec& x, rowvec& mean, mat& sigma, bool logd = false){ 
  
  // calculate density of multivariate normal distribution
  // args: x: row vector data
  //      mean: row vector mean, sigma: covariance matrix  
  //      logd: true for taking log
  // returns: out: pdf (or log pdf) of multivariate normal distribution
  
  int xdim = x.size(); 
  mat rooti = trans(inv(trimatu(chol(sigma))));
  double rootisum = sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0)*log2pi;
  
  vec z = rooti*trans(x-mean);
  double out = constants-0.5*sum(z%z)+rootisum;
  
  if (logd == false){ out = exp(out); }
  return(out);
}



// [[Rcpp::export]]
mat rmvn_rcpp(const int n, vec& mean, mat& sigma){
  
  // randomly generate samples from multivariate normal distribution
  // args: n: number of data 
  //      mean: row vector mean, sigma: covariance matrix  
  // returns: out: random samples from multivariate normal distribution
  
  int k = sigma.n_cols; // dimension of the multivariate normal distribution
  mat z = randn(n, k);
  mat out = repmat(mean,1,n).t()+z*chol(sigma);
  return(out);
}



// [[Rcpp::export]]
mat riwish_rcpp(const int df, mat& S){
  
  // randomly generate matrix from inverse wishart distribution
  // args: df: degrees of freedom
  //       S: inverse of scale matrix 
  // returns: out: random matrix from inverse wishart distribution
  
  S = S.i(); 
  int m = S.n_rows;
  
  mat Z(m, m); // Bartlett decomposition
  for (int i = 0; i < m; i++){
    Z(i,i) = sqrt(R::rchisq(df-i)); // Fill the diagonal
  }
  for (int j = 0; j < m; j++){  
    for(int i = j+1; i < m; i++){    
      Z(i,j) = R::rnorm(0,1); // Fill the lower matrix 
    }
  }
  
  mat C = trimatl(Z).t() * chol(S);
  // Random matrix from inverse wishart distribution
  mat out = (C.t()*C).i();
  return(out);
}



// [[Rcpp::export]]
mat kronecker_product_rcpp(mat& A, const int rowa, const int cola, mat& B, const int rowb, const int colb){
  
  // calculate the kronecker product of matrix A and B
  // args: A, B: matrix 
  //       rowa, cola, rowb, colb: number of rows and columns 
  // returns: C: A %x% B
  
  mat C(rowa*rowb, cola*colb);
  
  for (int i=0; i<rowa; i++){
    for (int k=0; k<rowb; k++){ 
      for (int j=0; j<cola; j++){ 
        for (int l=0; l<colb; l++){ 
          // Each element of matrix A is multiplied by 
          // whole Matrix B and stored as Matrix C 
          C(rowb*i+k,colb*j+l) = A(i,j) * B(k,l); 
        }
      }
    }
  }
  return(C);
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
double dmvn_kronecker_rcpp(const int xdim, vec& x, mat& sigmaA_inv, mat& sigmaB_inv, const int dimA, const int dimB){ 
  
  // calculate density of multivariate normal distribution (with kronecker product covariance matrix)
  // args: xdim: dimension of x, x: vector data
  //       sigmaA_inv, sigmaB_inv: inverse of cholesky decomposition of covariance matrix
  //       dimA, dimB: dimension of matrix A and B
  // returns: out: log pdf of multivariate normal distribution
  
  mat rooti = trans(kronecker_product_rcpp(sigmaA_inv, dimA, dimA, sigmaB_inv, dimB, dimB));
  double rootisum = sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0)*log2pi;
  vec z = rooti*x;
  double out = constants-0.5*sum(z%z)+rootisum;
  return(out);
}



// [[Rcpp::export]]
double rinvgamma_rcpp(double a, double b){
  
  // generate random samples from inverse-gamma distribution
  // args: inverse-gamma(a, b)
  // returns: random sample from inverse-gamma distribution
  
  return(1/R::rgamma(a, 1/b));
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



// [[Rcpp::export]]
mat update_theta_rcpp(const int n, const int M, const int S_tilde, const int Q, const int Nk_sum, 
                      vec& J, cube& X_tilde, cube& V, cube& U, cube& y, cube& alpha, mat& delta, 
                      cube& omega, vec& sigma2, mat& Sigma_theta_inv, mat& Sigma_theta_inv_mu_theta){
  
  // update theta
  // args: n, M, S_tilde, Q, Nk_sum, J, X_tilde, V, U, y: data
  //       alpha, delta, omega, sigma2: parameters
  //       Sigma_theta_inv, Sigma_theta_inv_mu_theta: hyper-parameters
  // returns: theta_update
  
  mat theta_update(M, S_tilde);
  vec mu_n(S_tilde); mat V_n(S_tilde, S_tilde);
  
  for (int m=0; m<M; m++){
    
    mat X_sum(S_tilde, S_tilde); X_sum.fill(0);
    vec Xy_sum(S_tilde); Xy_sum.fill(0);
    
    for (int i=0; i<n; i++){
      
      // pre-calculated values 
      mat X_tilde_i = X_tilde.subcube(i, 0, 0, i, J[i]-1, S_tilde-1);
      mat V_i = V.subcube(i, 0, 0, i, J[i]-1, Q-1);
      vec alpha_im = alpha.subcube(i, m, 0, i, m, Q-1);
      mat U_i = U.subcube(i, 0, 0, i, J[i]-1, Nk_sum-1);  
      vec delta_m = conv_to<vec>::from(delta.row(m));
      vec omega_im = omega.subcube(i, m, 0, i, m, J[i]-1);
      vec y_tilde_im = y.subcube(i, m, 0, i, m, J[i]-1);
      y_tilde_im -= V_i * alpha_im + U_i * delta_m + omega_im;
      
      // calculate X_sum and Xy_sum
      X_sum += X_tilde_i.t() * X_tilde_i;
      Xy_sum += X_tilde_i.t() * y_tilde_im;
    }
    
    // posterior mean and covariance matrix 
    V_n = inv_sympd(X_sum/sigma2[m] + Sigma_theta_inv);
    mu_n = V_n * (Xy_sum/sigma2[m] + Sigma_theta_inv_mu_theta);
    theta_update.row(m) = rmvn_rcpp(1, mu_n, V_n);
  }
  
  return(theta_update);
}



// [[Rcpp::export]]
cube update_alpha_rcpp(const int n, const int M, const int S_tilde, const int Q, const int Nk_sum, vec& J, cube& X_tilde, cube& V, 
                       cube& U, cube& y, mat& theta, mat& delta, cube& omega, vec& sigma2, cube& Sigma_alpha){
  
  // update alpha
  // args: n, M, S_tilde, Q, Nk_sum, J, X_tilde, V, U, y: data
  //       theta, delta, omega, sigma2, Sigma_alpha: parameters
  // returns: alpha_update
  
  cube alpha_update(n, M, Q);
  vec mu_n(Q); mat V_n(Q, Q);
  
  for (int m=0; m<M; m++){
    
    mat Sigma_alpha_m = Sigma_alpha.subcube(m, 0, 0, m, Q-1, Q-1);
    mat Sigma_alpha_m_inv = inv_sympd(Sigma_alpha_m);
    
    for (int i=0; i<n; i++){
      
      // pre-calculated values 
      mat V_i = V.subcube(i, 0, 0, i, J[i]-1, Q-1);
      mat X_tilde_i = X_tilde.subcube(i, 0, 0, i, J[i]-1, S_tilde-1);
      vec theta_m = conv_to<vec>::from(theta.row(m));
      mat U_i = U.subcube(i, 0, 0, i, J[i]-1, Nk_sum-1);  
      vec delta_m = conv_to<vec>::from(delta.row(m));
      vec omega_im = omega.subcube(i, m, 0, i, m, J[i]-1);
      vec y_tilde_im = y.subcube(i, m, 0, i, m, J[i]-1);
      y_tilde_im -= X_tilde_i * theta_m + U_i * delta_m + omega_im;
      
      // posterior mean and covariance matrix 
      V_n = inv_sympd(V_i.t() * V_i/sigma2[m] + Sigma_alpha_m_inv);
      mu_n = V_n * (V_i.t() * y_tilde_im/sigma2[m]);
      alpha_update.subcube(i,m,0,i,m,1) = rmvn_rcpp(1, mu_n, V_n);
    }
  }
  
  return(alpha_update);
}



// [[Rcpp::export]]
List update_delta_rcpp(const int n, const int m, const int S_tilde, const int Q, const int Nk_sum, vec& J, 
                       cube& X_tilde, cube& V, cube& U, cube& y, mat& theta, cube& alpha, cube& omega, vec& sigma2){
  
  // update posterior mean and covariance matrix for delta_m
  // args: n, m, S_tilde, Q, Nk_sum, J, X_tilde, V, U, y: data
  //       theta, alpha, omega, sigma2: parameters
  // returns: mu_n, V_n: posterior mean and covariance matrix for delta_m
  
  vec mu_n(Nk_sum); mat V_n(Nk_sum, Nk_sum);
  mat U_sum(Nk_sum, Nk_sum); U_sum.fill(0);
  vec Uy_sum(Nk_sum); Uy_sum.fill(0);
  
  for (int i=0; i<n; i++){
    
    // pre-calculated values 
    mat X_tilde_i = X_tilde.subcube(i, 0, 0, i, J[i]-1, S_tilde-1);
    vec theta_m = conv_to<vec>::from(theta.row(m));
    mat V_i = V.subcube(i, 0, 0, i, J[i]-1, Q-1);
    vec alpha_im = alpha.subcube(i, m, 0, i, m, Q-1);
    mat U_i = U.subcube(i, 0, 0, i, J[i]-1, Nk_sum-1);  
    vec omega_im = omega.subcube(i, m, 0, i, m, J[i]-1);
    vec y_tilde_im = y.subcube(i, m, 0, i, m, J[i]-1);
    y_tilde_im -= X_tilde_i * theta_m + V_i * alpha_im + omega_im;
    
    // calculate X_sum and Xy_sum
    U_sum += U_i.t() * U_i;
    Uy_sum += U_i.t() * y_tilde_im;
  }
  
  // posterior mean and covariance matrix 
  V_n = inv_sympd(U_sum/sigma2[m]);
  mu_n = V_n * (Uy_sum/sigma2[m]);
  
  // returnlist
  List returnlist;
  returnlist["mu_n"] = mu_n;
  returnlist["V_n"] = V_n;
  return(returnlist);
}



// [[Rcpp::export]]
cube update_omega_rcpp(const int n, const int M, const int S_tilde, const int Q, const int Nk_sum, 
                       vec& J, cube& X_tilde, cube& V, cube& U, cube& y, mat& t, 
                       mat& theta, mat delta, cube& alpha, mat& Cm, const double rho, vec& sigma2){
  
  // update omega
  // args: n, M, S_tilde, Q, Nk_sum, J, X_tilde, V, U, y, t: data
  //       theta, delta, alpha, Cm, rho, sigma2: parameters
  // returns: omega_update
  
  cube omega_update(n, M, max(J));
  
  mat D_inv = diagmat(1/sigma2); 
  mat Cm_inv = inv_sympd(Cm); 
  
  for (int i=0; i<n; i++){
    
    // calcualte the time-correlation matrix 
    vec t_i = conv_to<vec>::from(t.submat(i, 0, i, J[i]-1));
    mat Ct_i = auto_regressive_rcpp(J[i], 1, rho, t_i);
    
    // calculate the mean function
    int count = 0; vec mu_i(M*J[i]);
    for (int m=0; m<M; m++){
      for (int j=0; j<J[i]; j++){
        
        rowvec X_tilde_ij = X_tilde.subcube(i, j, 0, i, j, S_tilde-1);
        vec theta_m = conv_to<vec>::from(theta.row(m));
        rowvec V_ij = V.subcube(i, j, 0, i, j, Q-1);
        vec alpha_im = alpha.subcube(i, m, 0, i, m, Q-1);
        rowvec U_ij = U.subcube(i, j, 0, i, j, Nk_sum-1); 
        vec delta_m = conv_to<vec>::from(delta.row(m)); 
        
        double y_imj = y(i,m,j);
        double X_tilde_theta_imj = conv_to<double>::from(X_tilde_ij*theta_m);
        double V_alpha_imj = conv_to<double>::from(V_ij*alpha_im);
        double U_delta_imj = conv_to<double>::from(U_ij*delta_m);
        
        mu_i(count) = y_imj - X_tilde_theta_imj - V_alpha_imj - U_delta_imj;
        count += 1;
      }
    }
    
    // posterior mean and covariance matrix 
    mat I_i = 1.0 * eye<mat>(J[i], J[i]); mat Ct_i_inv = inv_sympd(Ct_i);
    mat Sigma_inv = kronecker_product_rcpp(D_inv, M, M, I_i, J[i], J[i]);
    mat Sigma_omega_inv = kronecker_product_rcpp(Cm_inv, M, M, Ct_i_inv, J[i], J[i]);
    mat V_n = inv_sympd(Sigma_omega_inv + Sigma_inv);
    vec mu_n = V_n * (Sigma_inv * mu_i);
    
    // update omgea_i, i = 1,2,...,n
    mat omega_i_update = rmvn_rcpp(1, mu_n, V_n);
    omega_i_update.reshape(J[i], M); omega_i_update = omega_i_update.t(); 
    omega_update.subcube(i,0,0,i,M-1,J[i]-1) = omega_i_update;
  }
  
  return(omega_update);
}



// [[Rcpp::export]]
double logll_Sigma_omega_rcpp(const int n, const int M, vec& J, mat & t, cube& omega, mat& Cm, const double rho){
  
  // calculate the log-likelihood of Sigma_omega
  // args: n, M, J, t: data
  //       omega, Cm, rho: parameters
  // returns: logll: log-likelihood of Sigma_omega
  
  double logll = 0;
  
  if (!(Cm.is_sympd())) { return(R_NegInf); }
  mat Cm_inv = inv(trimatu(chol(Cm)));
  
  for (int i=0; i<n; i++){
    
    // calcualte the time-correlation matrix 
    vec t_i = conv_to<vec>::from(t.submat(i,0,i,J[i]-1));
    mat Ct_i = auto_regressive_rcpp(J[i], 1, rho, t_i);
    mat Ct_i_inv = inv(trimatu(chol(Ct_i)));
    
    // calculate the loglikelihood of Sigma_omega_i, i = 1,2,...,n
    mat omega_i_mat = omega.subcube(i,0,0,i,M-1,J[i]-1);
    vec omega_i = omega_i_mat.t().as_col();
    logll += dmvn_kronecker_rcpp(M*J[i], omega_i, Cm_inv, Ct_i_inv, M, J[i]);
  }
  
  return(logll);
}



// [[Rcpp::export]]
double logll_Cm_rcpp(const int n, const int M, vec& J, cube& omega, mat& Am, vec& Bm, cube& Ct_inv){
  
  // calculate the log-likelihood of Cm
  // args: n, M, J: data
  //       omega, Am, Bm, Ct_inv: parameters
  // returns: logll: log-likelihood of Cm
  
  double logll = 0;
  
  mat Cm = Am * Am.t() + diagmat(Bm);
  if (!(Cm.is_sympd())) { return(R_NegInf); }
  mat Cm_inv = inv(trimatu(chol(Cm)));
  
  for (int i=0; i<n; i++){
    
    // calculate the loglikelihood of Sigma_omega_i, i = 1,2,...,n
    mat Ct_i_inv = Ct_inv.subcube(i,0,0,i,J[i]-1,J[i]-1);
    mat omega_i_mat = omega.subcube(i,0,0,i,M-1,J[i]-1);
    vec omega_i = omega_i_mat.t().as_col();
    logll += dmvn_kronecker_rcpp(M*J[i], omega_i, Cm_inv, Ct_i_inv, M, J[i]);
  }
  
  return(logll);
}



// [[Rcpp::export]]
List logll_rho_rcpp(const int n, const int M, vec& J, mat & t, cube& omega, mat& Cm, const double rho){
  
  // calculate the log-likelihood of Sigma_omega
  // args: n, M, J, t: data
  //       omega, Cm, rho: parameters
  // returns: logll: log-likelihood of Sigma_omega
  
  double logll = 0; 
  cube Ct_inv(n, max(J), max(J));
  
  if (!(Cm.is_sympd())) { return(R_NegInf); }
  mat Cm_inv = inv(trimatu(chol(Cm)));
  
  for (int i=0; i<n; i++){
    
    // calcualte the time-correlation matrix 
    vec t_i = conv_to<vec>::from(t.submat(i,0,i,J[i]-1));
    mat Ct_i = auto_regressive_rcpp(J[i], 1, rho, t_i);
    mat Ct_i_inv = inv(trimatu(chol(Ct_i)));
    Ct_inv.subcube(i,0,0,i,J[i]-1,J[i]-1) = Ct_i_inv;
    
    // calculate the loglikelihood of Sigma_omega_i, i = 1,2,...,n
    mat omega_i_mat = omega.subcube(i,0,0,i,M-1,J[i]-1);
    vec omega_i = omega_i_mat.t().as_col();
    logll += dmvn_kronecker_rcpp(M*J[i], omega_i, Cm_inv, Ct_i_inv, M, J[i]);
  }
  
  // returnlist
  List returnlist;
  returnlist["logll"] = logll;
  returnlist["Ct_inv"] = Ct_inv;
  return(returnlist);
}



// [[Rcpp::export]]
vec update_sigma2_rcpp(const int n, const int M, const int S_tilde, const int Q, const int Nk_sum, vec& J, cube& X_tilde, cube& V, 
                       cube& U, cube& y, mat& theta, mat& delta, cube& alpha, cube& omega, const double g_1, const double g_2){
  
  // update sigma2
  // args: n, M, S_tilde, Q, Nk_sum, J, X_tilde, V, U, y: data
  //       theta, delta, alpha, omega: parameters
  //       g_1, g_2: hyper-parameters
  // returns: sigma2_update
  
  vec sigma2_udpate(M);
  
  for (int m=0; m<M; m++){
    
    double y_tilde_sum = 0;
    
    for (int i=0; i<n; i++){
      
      // pre-calculated values 
      mat X_tilde_i = X_tilde.subcube(i, 0, 0, i, J[i]-1, S_tilde-1);
      vec theta_m = conv_to<vec>::from(theta.row(m));
      mat V_i = V.subcube(i, 0, 0, i, J[i]-1, Q-1);
      vec alpha_im = alpha.subcube(i, m, 0, i, m, Q-1);
      mat U_i = U.subcube(i, 0, 0, i, J[i]-1, Nk_sum-1);  
      vec delta_m = conv_to<vec>::from(delta.row(m));
      vec omega_im = omega.subcube(i, m, 0, i, m, J[i]-1);
      vec y_tilde_im = y.subcube(i, m, 0, i, m, J[i]-1);
      y_tilde_im -= X_tilde_i * theta_m + V_i * alpha_im + U_i * delta_m + omega_im;
      
      // calculate y_tilde_sum
      for (int j=0; j<J[i]; j++){ y_tilde_sum += pow(y_tilde_im(j), 2); }
    }
    
    // inverse-gamma posterior distribution
    double g_1_star = g_1 + accu(J)/2;
    double g_2_star = g_2 + y_tilde_sum/2;
    sigma2_udpate[m] = rinvgamma_rcpp(g_1_star, g_2_star);
  }
  
  return(sigma2_udpate);
}



// [[Rcpp::export]]
cube update_y_rcpp(const int n, const int M, const int S_tilde, const int Q, const int Nk_sum,
                   cube& data_index, vec& J, cube& X_tilde, cube& V, cube& U, cube& y, mat& t, 
                   mat& theta, mat& delta, cube& alpha, mat& Cm, const double rho, vec& sigma2){
  
  // update y
  // args: n, M, S_tilde, Q, Nk_sum, data_index, J, X_tilde, V, U, y, t: data
  //       theta, delta, alpha, Cm, rho, sigma2: parameters
  // returns: y_update
  
  cube y_update = y;
  mat D = diagmat(sigma2); 
  
  for (int i=0; i<n; i++){
    
    // calcualte the time-correlation matrix 
    vec t_i = conv_to<vec>::from(t.submat(i,0,i,J[i]-1));
    mat Ct_i = auto_regressive_rcpp(J[i], 1, rho, t_i);
    
    // calculate full covariance matrix 
    mat I_i = 1.0 * eye<mat>(J[i], J[i]); 
    mat Sigma_i = kronecker_product_rcpp(D, M, M, I_i, J[i], J[i]);
    Sigma_i += kronecker_product_rcpp(Cm, M, M, Ct_i, J[i], J[i]);
    
    // index of missing and non-missing data 
    mat data_index_i = data_index.subcube(i, 0, 0, i, M-1, J[i]-1);
    uvec update_index = find(data_index_i.t() == 0);
    uvec obs_index = find(data_index_i.t() == 1);
    
    // calculate the mean function
    int count1 = 0; vec mu_i(M*J[i]);
    for (int m=0; m<M; m++){
      for (int j=0; j<J[i]; j++){
          
        rowvec X_tilde_ij = X_tilde.subcube(i, j, 0, i, j, S_tilde-1);
        vec theta_m = conv_to<vec>::from(theta.row(m));
        rowvec V_ij = V.subcube(i, j, 0, i, j, Q-1);
        vec alpha_im = alpha.subcube(i, m, 0, i, m, Q-1);
        rowvec U_ij = U.subcube(i, j, 0, i, j, Nk_sum-1); 
        vec delta_m = conv_to<vec>::from(delta.row(m)); 
          
        double X_tilde_theta_imj = conv_to<double>::from(X_tilde_ij*theta_m);
        double V_alpha_imj = conv_to<double>::from(V_ij*alpha_im);
        double U_delta_imj = conv_to<double>::from(U_ij*delta_m);
        
        if (data_index_i(m,j) == 1){  
          double y_imj = y_update(i,m,j);
          mu_i(count1) = y_imj - X_tilde_theta_imj - V_alpha_imj - U_delta_imj;
        }else{ mu_i(count1) = X_tilde_theta_imj + V_alpha_imj + U_delta_imj; }
        count1 += 1;
      }
    }
    vec mu_x = mu_i.elem(update_index);
    vec mu_y = mu_i.elem(obs_index);
    
    // calculate covariance matrix 
    mat Sigma_xx = Sigma_i.submat(update_index,update_index);
    mat Sigma_xy = Sigma_i.submat(update_index,obs_index);
    mat Sigma_yx = Sigma_i.submat(obs_index,update_index);
    mat Sigma_yy = Sigma_i.submat(obs_index,obs_index);
    
    // posterior mean and covariance matrix 
    vec mu = mu_x + Sigma_xy * inv_sympd(Sigma_yy) * mu_y;
    mat Sigma = Sigma_xx - Sigma_xy * inv_sympd(Sigma_yy) * Sigma_yx;
    
    // update omgea_i, i = 1,2,...,n
    int count2 = 0;
    vec y_i_update = conv_to<vec>::from(rmvn_rcpp(1, mu, Sigma));
    for (int m=0; m<M; m++){
      for (int j=0; j<J[i]; j++){
        if (data_index_i(m,j) == 0){
          y_update(i,m,j) = y_i_update(count2);
          count2 += 1;
        }
      }
    }
  }
  
  return(y_update);
}