####################################################### MCMC ########################################################


MCMC_Main <- function(Nit, burn.in, thin.fac){
  
  # MCMC main function 
  # args: Nit: number of intertaions 
  #       burn.in: the burn-in iterations
  #       thin.fac: thinning factor for post burn-in samples 
  # returns: MCMC posterior samples after burn-in with thinning factor 
  
  mcmc <- NULL
  mcmc$y <- array(NA, dim=c(Nit, n, M, max(J)))
  mcmc$theta <- array(NA, dim=c(Nit, M, S_tilde))
  mcmc$delta <- array(NA, dim=c(Nit, M, Nk_sum))
  mcmc$alpha <- array(NA, dim=c(Nit, n, M, Q))
  mcmc$Sigma_alpha <- array(NA, dim=c(Nit, M, Q, Q))
  mcmc$Am <- array(NA, dim=c(Nit, M, L))
  mcmc$Bm <- matrix(NA, nrow=Nit, ncol=M)
  mcmc$Cm <- array(NA, dim=c(Nit, M, M))
  mcmc$rho <- rep(NA, Nit)
  mcmc$sigma2 <- matrix(NA, nrow=Nit, ncol=M)
  mcmc$lambda <- matrix(NA, nrow=Nit, ncol=S_minus)
  mcmc$tau <- rep(NA, Nit)
  mcmc$nu <- matrix(NA, nrow=Nit, ncol=S_minus)
  mcmc$psi <- rep(NA, Nit)
  
  ### Initialize
  initial <- init()
  mcmc$y[1,,,] <- y
  mcmc$theta[1,,] <- initial$theta
  mcmc$delta[1,,] <- initial$delta
  mcmc$alpha[1,,,] <- initial$alpha
  mcmc$Sigma_alpha[1,,,] <- initial$Sigma_alpha
  mcmc$Am[1,,] <- initial$Am
  mcmc$Bm[1,] <- initial$Bm
  mcmc$Cm[1,,] <- initial$Cm
  mcmc$rho[1] <- initial$rho
  mcmc$sigma2[1,] <- initial$sigma2
  mcmc$omega <- initial$omega
  mcmc$Ct_inv <- initial$Ct_inv
  mcmc$lambda[1,] <- initial$lambda
  mcmc$tau[1] <- initial$tau
  mcmc$nu[1,] <- initial$nu
  mcmc$psi[1] <- initial$psi
  mcmc$Sigma_theta_inv <- initial$Sigma_theta_inv
  mcmc$Sigma_theta_inv_mu_theta <- initial$Sigma_theta_inv_mu_theta
  
  ### Start of the chain
  
  start.time = proc.time()
  
  for (nit in 2:Nit){
    
    print(nit)
    
    # update y
    mcmc$y[nit,,,] <- update_y_rcpp(n, M, S_tilde, Q, Nk_sum, data_index, J, X_tilde, V, U, mcmc$y[nit-1,,,], t, mcmc$theta[nit-1,,], mcmc$delta[nit-1,,], mcmc$alpha[nit-1,,,], mcmc$Cm[nit-1,,], mcmc$rho[nit-1], mcmc$sigma2[nit-1,])
    
    # update theta
    mcmc$theta[nit,,] <- update_theta_rcpp(n, M, S_tilde, Q, Nk_sum, J, X_tilde, V, U, mcmc$y[nit,,,], mcmc$alpha[nit-1,,,], mcmc$delta[nit-1,,], mcmc$omega, mcmc$sigma2[nit-1,], mcmc$Sigma_theta_inv, mcmc$Sigma_theta_inv_mu_theta)
    
    # update lambda, tau, nu, psi, Sigma_theta_inv
    mcmc$lambda[nit,] <- update_lambda(mcmc$theta[nit,,shrink_index], mcmc$tau[nit-1], mcmc$nu[nit-1,])
    mcmc$tau[nit] <- update_tau(mcmc$theta[nit,,shrink_index], mcmc$lambda[nit,], mcmc$psi[nit-1])
    mcmc$nu[nit,] <- update_nu(mcmc$lambda[nit,])
    mcmc$psi[nit] <- update_psi(mcmc$tau[nit])
    mcmc$Sigma_theta_inv[shrink_index,shrink_index] <- 1/(mcmc$tau[nit]^2)*diag(1/mcmc$lambda[nit,]^2)
    
    # update alpha
    mcmc$alpha[nit,,,] <- update_alpha_rcpp(n, M, S_tilde, Q, Nk_sum, J, X_tilde, V, U, mcmc$y[nit,,,], mcmc$theta[nit,,], mcmc$delta[nit-1,,], mcmc$omega, mcmc$sigma2[nit-1,], mcmc$Sigma_alpha[nit-1,,,])
    
    # update Sigma_alpha
    mcmc$Sigma_alpha[nit,,,] <- update_Sigma_alpha(mcmc$alpha[nit,,,])
    
    # update delta
    mcmc$delta[nit,,] <- update_delta(mcmc$theta[nit,,], mcmc$alpha[nit,,,], mcmc$omega, mcmc$sigma2[nit-1,], mcmc$y[nit,,,])
    
    # update omega
    mcmc$omega <- update_omega_rcpp(n, M, S_tilde, Q, Nk_sum, J, X_tilde, V, U, mcmc$y[nit,,,], t, mcmc$theta[nit,,], mcmc$delta[nit,,], mcmc$alpha[nit,,,], mcmc$Cm[nit-1,,], mcmc$rho[nit-1], mcmc$sigma2[nit-1,])
    
    # update Cm
    mcmc$Am[nit,,] <- update_Am(mcmc$omega, mcmc$Am[nit-1,,], mcmc$Bm[nit-1,], mcmc$Ct_inv)
    mcmc$Bm[nit,] <- update_Bm(mcmc$omega, mcmc$Am[nit,,], mcmc$Bm[nit-1,], mcmc$Ct_inv)
    mcmc$Cm[nit,,] <- mcmc$Am[nit,,] %*% t(mcmc$Am[nit,,]) + diag(mcmc$Bm[nit,])
    
    # update rho and Ct_inv
    rho_update_list <- update_rho(mcmc$omega, mcmc$Cm[nit,,], mcmc$rho[nit-1], mcmc$Ct_inv)
    mcmc$rho[nit] <- rho_update_list$rho; mcmc$Ct_inv <- rho_update_list$Ct_inv
    
    # update sigma2
    mcmc$sigma2[nit,] <- update_sigma2_rcpp(n, M, S_tilde, Q, Nk_sum, J, X_tilde, V, U, mcmc$y[nit,,,], mcmc$theta[nit,,], mcmc$delta[nit,,], mcmc$alpha[nit,,,], mcmc$omega, g_1, g_2)                        
  }
  
  duration = proc.time()-start.time
  print(duration) # print the running time 
  
  ### End of the chain
  
  ### Posterior samples
  post_index <- seq(burn.in+1, Nit, by=thin.fac) # index of posterior samples
  
  post <- NULL
  post$y <- mcmc$y[post_index,,,]
  post$theta <- mcmc$theta[post_index,,]
  post$delta <- mcmc$delta[post_index,,]
  post$alpha <- mcmc$alpha[post_index,,,]
  post$Sigma_alpha <- mcmc$Sigma_alpha[post_index,,,]
  post$Am <- mcmc$Am[post_index,,]
  post$Bm <- mcmc$Bm[post_index,]
  post$Cm <- mcmc$Cm[post_index,,]
  post$rho <- mcmc$rho[post_index]
  post$sigma2 <- mcmc$sigma2[post_index,]
  post$lambda <- mcmc$lambda[post_index,]
  post$tau <- mcmc$tau[post_index]
  post$chi <- 1/(1+mean(post$tau)^2*colMeans(post$lambda)^2)
  
  return(post)
}
