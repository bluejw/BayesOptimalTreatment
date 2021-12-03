################################################## MCMC Functions ##############################################



init <- function(){
  
  # MCMC initialization 
  # initialize parameters need to be estimated in MCMC
  
  theta_init <- rmvn_rcpp(M, rep(0,S_tilde), diag(1,S_tilde))
  Sigma_alpha_init <- array(NA, dim=c(M, Q, Q))
  delta_init <- matrix(0.5, nrow=M, ncol=Nk_sum)
  for (m in 1:M) { Sigma_alpha_init[m,,] <- riwish_rcpp(Q+1, diag(1,Q)) }
  alpha_init <- array(NA, dim=c(n, M, Q))
  for (i in 1:n){ for (m in 1:M){ alpha_init[i,m,] <- rmvn_rcpp(1, rep(0,Q), Sigma_alpha_init[m,,]) }}
  
  Am_init <- matrix(1, nrow=M, ncol=L)
  Bm_init <- rep(1, M)
  Cm_init <- Am_init %*% t(Am_init) + diag(Bm_init)
  rho_init <- 0.5 
  sigma2_init <- rep(1, M)
  omega_init <- array(NA, dim=c(n, M, max(J))) 
  Ct_inv_init <- array(NA, dim=c(n, max(J), max(J)))
  for (i in 1:n){
    Ct_i_init <- matrix(NA, nrow=J[i], ncol=J[i])
    for (j in 1:J[i]){ Ct_i_init[j,1:J[i]] <- rho_init^(abs(t[i,j]-t[i,1:J[i]])) }
    Ct_inv_init[i,1:J[i],1:J[i]] <- solve(chol(Ct_i_init))
    omega_i <- rmvn_rcpp(1, rep(0,M*J[i]), Cm_init %x% Ct_i_init)
    omega_init[i,1:M,1:J[i]] <- matrix(omega_i, nrow=M, ncol=J[i], byrow=T)
  }
  
  lambda_init <- rep(1, S_minus); tau_init <- 1
  nu_init <- rep(1, S_minus); psi_init <- 1
  Sigma_theta_inv_init <- diag(1, S_tilde)
  Sigma_theta_inv_init[1:(S+D),1:(S+D)] <- diag(1/100, S+D); 
  Sigma_theta_inv_init[shrink_index,shrink_index] <- 1/(tau_init^2)*diag(1/lambda_init^2)
  mu_theta_init <- rep(0, S_tilde); Sigma_theta_inv_mu_theta_init <- Sigma_theta_inv_init %*% mu_theta_init
  
  init_list <- list(theta=theta_init, delta=delta_init, alpha=alpha_init,Sigma_alpha=Sigma_alpha_init, omega=omega_init, 
                    Cm=Cm_init, rho=rho_init, Am=Am_init, Bm=Bm_init, Cm=Cm_init, Ct_inv=Ct_inv_init, sigma2=sigma2_init,
                    lambda=lambda_init, tau=tau_init, nu=nu_init, psi=psi_init, 
                    Sigma_theta_inv=Sigma_theta_inv_init, Sigma_theta_inv_mu_theta=Sigma_theta_inv_mu_theta_init)
  return(init_list)
}



update_lambda <- function(theta, tau, nu){
  
  # update lambda
  # args: theta: M*S_minus dimensional parameters with shrinkage priors
  #       tau, nu: hyper-parameters
  # returns: lambda_update
  
  lambda_update <- rep(NA, S_minus)
  
  for (s in 1:S_minus){
    a_star <- (M+1)/2
    b_star <- 1/nu[s] + sum(theta[,s]^2)/(2*tau^2)
    lambda_update[s] <- sqrt(rinvgamma(1, a_star, b_star))
  }
  return(lambda_update)
}



update_tau <- function(theta, lambda, psi){
  
  # update tau
  # args: theta: M*S_minus dimensional parameters with shrinkage priors
  #       lambda, psi: hyper-parameters
  # returns: tau_update
  
  theta_star <- 0
  for (s in 1:S_minus){ 
    theta_star <- theta_star + sum((theta[,s]/lambda[s])^2) 
  }
  a_star <- (M*S_minus+1)/2
  b_star <- 1/psi + theta_star/2
  tau_update <- sqrt(rinvgamma(1, a_star, b_star))
  return(tau_update)
}



update_nu <- function(lambda){
  
  # update nu
  # args: lambda: hyper-parameters
  # returns: nu_update
  
  nu_update <- rep(NA, S_minus)
  for (s in 1:S_minus){
    nu_update[s] <- rinvgamma(1, 1, 1+1/lambda[s]^2)
  }
  return(nu_update)
}



update_psi <- function(tau){
  
  # update psi
  # args: tau: hyper-parameters
  # returns: psi_update
  
  psi_update <- rinvgamma(1, 1+1/tau^2)
  return(psi_update)
}



update_Sigma_alpha <- function(alpha){
  
  # update Sigma_alpha
  # args: alpha: parameter
  # returns: Sigma_alpha_update
  
  Sigma_alpha_update <- array(NA, dim=c(M, Q, Q))
  
  for (m in 1:M){
    alpha_sum <- matrix(0, nrow=Q, ncol=Q)
    for (i in 1:n){ 
      alpha_sum <- alpha_sum + alpha[i,m,] %*% t(alpha[i,m,])
    }
    # Inverse Wishart posterior distribution 
    a_0_star <- a_0 + n; A_0_star <- A_0 + alpha_sum
    Sigma_alpha_update[m,,] <- riwish_rcpp(a_0_star, A_0_star) 
  }
  
  return(Sigma_alpha_update)
}



update_delta <- function(theta, alpha, omega, sigma2, y){
  
  # update delta
  # args: theta, alpha, omega, sigma2: parameters 
  # returns: delta_update
  
  delta_update <- matrix(NA, nrow=M, ncol=Nk_sum)
  
  for (m in 1:M){
    
    # calculate U_sum and Uy_sum
    update_list <- update_delta_rcpp(n, m-1, S_tilde, Q, Nk_sum, J, X_tilde, V, U, y, theta, alpha, omega, sigma2)
    mu_n <- update_list$mu_n; V_n <- update_list$V_n
    
    # multivariate truncated normal posterior distribution
    if (m %in% negative_state){ delta_update[m,] <- rtmvnorm(n=1, mu=mu_n, sigma=V_n, lb=rep(-Inf,Nk_sum), ub=rep(0,Nk_sum))
    }else{ delta_update[m,] <- rtmvnorm(n=1, mu=mu_n, sigma=V_n, lb=rep(0,Nk_sum), ub=rep(Inf,Nk_sum)) }
  }
  return(delta_update)
}



update_Am <- function(omega, Am, Bm, Ct_inv){
  
  # update Am
  # args: omega, Am, Bm, Ct_inv: parameters
  # returns: Am_update
  
  if (L == 1) { Am_update <- matrix(Am) 
  }else { Am_update <- Am }
  step <- 0.05 # step size of proposal distribution
  
  for (i in 1:M){
    for (j in 1:L){
      # lower and upper bounds for the (i,j) parameter
      lower <- Am_update[i,j]-step
      upper <- Am_update[i,j]+step
      # propose the new (i,j) parameter 
      Am_new <- Am_update
      Am_new[i,j] <- runif(1, lower, upper)
      # acceptance ratio with uniform prior and proposal 
      ratio <- logll_Cm_rcpp(n, M, J, omega, Am_new, Bm, Ct_inv) - logll_Cm_rcpp(n, M, J, omega, Am_update, Bm, Ct_inv)
      # accept or reject 
      if (log(runif(1)) < ratio){ Am_update <- Am_new }
    }
  }
  
  return(Am_update)
}



update_Bm <- function(omega, Am, Bm, Ct_inv){
  
  # update Bm
  # args: omega, Am, Bm, Ct_inv: parameters
  # returns: Bm_update
  
  Bm_update <- Bm
  if (L == 1) { Am <- matrix(Am) }
  step <- 0.05 # step size of proposal distribution
  
  for (i in 1:M){
    # lower and upper bounds for the i-th parameter
    lower <- Bm_update[i]-step
    upper <- Bm_update[i]+step
    # propose the new i-th parameter 
    Bm_new <- Bm_update
    Bm_new[i] <- runif(1, lower, upper)
    # acceptance ratio with uniform prior and proposal 
    ratio <- logll_Cm_rcpp(n, M, J, omega, Am, Bm_new, Ct_inv) - logll_Cm_rcpp(n, M, J, omega, Am, Bm_update, Ct_inv)
    # accept or reject 
    if (log(runif(1)) < ratio){ Bm_update <- Bm_new }
  }
  
  return(Bm_update)
}



update_rho <- function(omega, Cm, rho, Ct_inv){
  
  # update rho
  # args: omega, sigma2_omega, Cm, rho: parameters
  # returns: rho_update
  
  rho_update <- rho
  Ct_inv_update <- Ct_inv
  step <- 0.02 # step size of proposal distribution
  eps <- 1e-300 # avoid numerical issue
  
  # lower and upper bounds for the parameter
  lower <- max(rho_update-step, eps)
  upper <- min(rho_update+step, 1-eps)
  # propose the new parameter 
  rho_new <- runif(1, lower, upper)
  # acceptance ratio with uniform prior and proposal 
  update_list <- logll_rho_rcpp(n, M, J, t, omega, Cm, rho_new)
  ratio <- update_list$logll - logll_Sigma_omega_rcpp(n, M, J, t, omega, Cm, rho_update)
  # accept or reject 
  if (log(runif(1)) < ratio){ rho_update <- rho_new; Ct_inv_update <- update_list$Ct_inv }
  
  # returnlist 
  returnlist <- list(rho=rho_update, Ct_inv=Ct_inv_update)
  return(returnlist)
}