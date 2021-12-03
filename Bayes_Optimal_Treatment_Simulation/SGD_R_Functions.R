################################################## SGD Functions ##############################################



softmax1 <- function(k, y, theta1k){
  
  # calculate the probability for number of single drugs ck=0,1,..,Ck in drug class k
  # args: k: drug class index, y: state data, theta1k: level1 regimen model parameter 
  # returns: softmax probability
  
  temp <- rep(NA, Ck[k]+1); temp[1] <- 0
  if (Ck[k] > 0){
    for (ck in 2:(Ck[k]+1)){
      temp[ck] <- t(y) %*% theta1k[[ck-1]]
    }
  }
  
  exptemp <- exp(temp-max(temp))
  return(exptemp/sum(exptemp))
}



softmax2 <- function(k, y, theta2k){
  
  # calculate the probability for single drug usage in drug class k
  # args: theta2k: level2 regimen model parameter
  # returns: softmax probability
  
  temp <- rep(NA, Nk[k])
  for (nk in 1:Nk[k]){
    temp[nk] <- t(y) %*% theta2k[[nk]]
  }

  exptemp <- exp(temp-max(temp))
  return(exptemp/sum(exptemp))
}



regimen_tree <- function(Z){
  
  # represent the ART regimen in subset-tree structure 
  # args: Z: a vector of ART regimen, e.g., c("3TC", "AZT", "NFV")
  # returns: Ztree: a list of subset-tree representation of ART regimen Z
  
  # level 1
  a1 <- NULL # number of single drugs used in each drug class 
  a1$NRTI1 <- sum(Z %in% all_drugs$NRTI1)
  a1$NRTI2 <- sum(Z %in% all_drugs$NRTI2)
  a1$NNRTI <- sum(Z %in% all_drugs$NNRTI)
  a1$PI <- sum(Z %in% all_drugs$PI)
  a1$INSTI <- sum(Z %in% all_drugs$INSTI)
  a1$EI <- sum(Z %in% all_drugs$EI)
  
  # level 2
  a2 <- NULL # indicator of all possible single drug usuage in each drug class
  a2$NRTI1 <- as.numeric(all_drugs$NRTI1 %in% Z) 
  a2$NRTI2 <- as.numeric(all_drugs$NRTI2 %in% Z) 
  a2$NNRTI <- as.numeric(all_drugs$NNRTI %in% Z) 
  a2$PI <- as.numeric(all_drugs$PI %in% Z) 
  a2$INSTI <- as.numeric(all_drugs$INSTI %in% Z) 
  a2$EI <- as.numeric(all_drugs$EI %in% Z) 
  
  # ART regimen Z in subset-tree structure
  Ztree <- NULL
  Ztree$a1 <- a1
  Ztree$a2 <- a2
  return(Ztree)
}



grad_logll_regimen <- function(Zn, y, theta){
  
  # calculate the gradient of log-likelihood for the new regimen Zn 
  # args: Zn: regimens at j+1-th visit
  #       y: state matrix at j-th visit
  #       theta: regimen model parameters 
  # returns: grad
  
  # gradients
  grad <- NULL
  grad$theta1 <- NULL
  grad$theta2 <- NULL
  
  # parameters
  theta1 <- theta$theta1
  theta2 <- theta$theta2
  
  # add intercept to state value
  y <- c(1, (y-y_mean)/y_std)
  
  # reconstruct ART regimen in subset-tree structure
  Zntree <- regimen_tree(Zn)
  a1 <- Zntree$a1; a2 <- Zntree$a2
  
  # calculate the gradient
  for (k in 1:K){
    
    # gradients of theta1
    grad$theta1[[k]] <- list()
    a1k <- a1[[k]] # number of single drugs used in drug class k
    p1k <- softmax1(k, y, theta1[[k]]) # probability for number of single drugs in each drug class
    for (ck in 1:Ck[k]){
      if (ck == a1k){ grad$theta1[[k]][[ck]] <- (1-p1k[ck+1])*y 
      }else{ grad$theta1[[k]][[ck]] <- (-p1k[ck+1])*y }
    }
    
    # gradients of theta2
    if (k == K){
      # if it is EI drug class (or any drug class contains only one possible drug)
      grad$theta2[[k]] <- list(rep(0, M_tilde)) 
    }else{
      grad$theta2[[k]] <- list()
      a2k <- a2[[k]] # single drug indicator in drug class k
      omegak <- softmax2(k, y, theta2[[k]]) # odds/weights for single drug usage in each class
      p2k <- dMWNCHypergeo(x=a2k, m=rep(1,Nk[k]), n=a1k, odds=omegak) # probability for single drug usage in each class
      tauk <- omegak/sum(omegak*(1-a2k)) # Nk-dimensional vector in the density integrand 
        
      for (nk in 1:Nk[k]){
        # pre-calculated values for integral
        omegak_dev <- -omegak[nk]*omegak # derivative of omegak on theta2k[nk]
        omegak_dev[nk] <- omegak_dev[nk]+omegak[nk]
        omegak_sum <- sum(omegak*(1-a2k))
        omegak_devsum <- sum(omegak_dev*(1-a2k))
        grad$theta2[[k]][[nk]] <- rep(0, M_tilde)
        # calculate the gradient for theta2k[nk]
        for (m in 1:M_tilde){
          grad$theta2[[k]][[nk]][m] <- 1/p2k * 
            noncentral_hypergeo_integrate_rcpp(num=10, lower=0, upper=1,
            Nk[k], a2k, tauk, omegak, omegak_dev, omegak_sum, omegak_devsum, y[m])
        }
      }
    }
  }
  
  return(grad)
}



generate_regimen <- function(y, theta){
  
  # generate the new regimen Zn conditional on state y 
  # args: y: state matrix at j-th visit
  #       theta: regimen model parameters 
  # returns: Zn: regimen at j+1-th visits
  
  Zn <- c()
  
  # parameters
  theta1 <- theta$theta1
  theta2 <- theta$theta2

  # add intercept to state value
  y <- c(1, (y-y_mean)/y_std)
  
  # generate new regimen
  for (k in 1:K){
    p1k <- softmax1(k, y, theta1[[k]]) # probability for number of single drugs in each class
    if (k == 1 | k == 2){
      # if it is NRTI drug class
      a1k <- 1
    }else{
      # if it is not NRTI drug class
      a1k <- which(rmultinom(1, 1, prob=p1k)==1) - 1
    }
    omegak <- softmax2(k, y, theta2[[k]]) # odds/weights for single drug usage in each class
    a2k <- rMWNCHypergeo(nran=1, m=rep(1,Nk[k]), n=a1k, odds=omegak)
    Zn <- c(Zn, all_drugs[[k]][which(a2k==1)])
    if (k == 4 & a1k > 0) { Zn <- c(Zn, "RTV") } # include RTV when PI drugs selected
    if (k == 5 & c("EVG") %in% Zn) { Zn <- c(Zn, "COBI") } # include COBI when INSTI drug EVG selected
  }
  Zn <- sort(Zn)
  
  return(Zn)
}



generate_state <- function(data, jp, Zn, tnt, phi_nit, theta_all, delta_all, alpha_i_all){
  
  # generate the new state yn 
  # args: data: list of state, covarites, drug usuage, time at 1:jp-th visits
  #       Zn, tnt: regimen and time at jp+1-th visit
  #       phi_nit: posterior samples of state model parameters at one iteration
  # returns: yn: state matrix at jp+1-th visit
  #          data: the updated data matrix
  #          uq: uncertainty quantification
  
  yn <- rep(NA, M)
  
  # data
  yp <- data$yp
  Xp <- data$Xp
  Hp <- data$Hp
  XHp <- data$XHp
  Vp <- data$Vp
  tp <- data$tp
  Zp <- data$Zp
  Up <- data$Up
  Dp <- data$Dp
  
  # parameters
  theta <- phi_nit$theta
  delta <- phi_nit$delta
  alpha_i <- phi_nit$alpha_i
  Cm <- phi_nit$Cm
  rho <- phi_nit$rho
  sigma2 <- phi_nit$sigma2
  
  # time correlation matrix
  t_i <- c(tp, tnt) # time points for i-th subject
  Ct_i <- auto_regressive_rcpp(n=jp+1, sigma2=1, rho, t_i)
  
  # drug toxicity
  Z_i <- rbind(Zp, c(Zn, rep(NA, max_num_drugs-length(Zn))))
  D_i <- rbind(Dp, as.numeric(drug_names %in% Zn))
  tit <- t_i - t_i[1] + 0.5; tmax <- max(tit)
  drug_toxic <- as.vector(drug_toxic_integrate_rcpp(num=100, lower=0, upper=tmax, Nk_sum, jp+1, D_i, tit))
  U_i <- rbind(Up, (drug_toxic-U_mean)/U_std)
  
  # fixed and random effect covariates
  X_i <- rbind(Xp, Xp[1,1:S]) 
  V_i <- rbind(Vp, c(1, (tnt-time_mean)/time_std))
  drug_simi <- Drug_Similarity(Zn, z[index_kernel,]) # drug similarity between Zn and kernel knots
  if (pca){
    # if pca on kernel regression 
    if (sum(drug_simi) != 0){ 
      drug_simi <- drug_simi/sum(drug_simi)
      drug_simi <- (drug_simi-H_mean)/H_std
      drug_simi <- drug_simi %*% eigenvec[,1:D] 
      H_i <- rbind(Hp, (drug_simi-Hpca_mean)/Hpca_std)
    }else{ H_i <- rbind(Hp, drug_simi[1:D]) }
  }else{
    # if no pca on kernel regression
    if (sum(drug_simi) != 0){ 
      drug_simi <- drug_simi/sum(drug_simi)
      H_i <- rbind(Hp, (drug_simi-H_mean)/H_std)
    }else{ H_i <- rbind(Hp, drug_simi) }
  }
  
  # drug comb and baseline interaction terms
  drug_base_inter <- Xp[1,2:S] %x% as.vector(drug_simi)
  XH_i <- rbind(XHp, (drug_base_inter-XH_mean)/XH_std)
  X_tilde_i <- cbind(X_i, H_i, XH_i)
  
  # mean and covariance matrix for Gaussian process
  mu_i <- matrix(NA, nrow=M, ncol=jp+1)
  Sigma_i <- Cm[1:M,1:M] %x% Ct_i + diag(sigma2, M) %x% diag(1, jp+1) 
  mu_i_all <- matrix(NA, nrow=M, ncol=num_samples)
  for (m in 1:M){
    # generate the new state matrix
    mu_i[m,] <- X_tilde_i %*% theta[m,] + U_i %*% delta[m,] + V_i %*% alpha_i[m,]
    start <- (m-1)*(jp+1)+1; end <- m*(jp+1); index <- 1:jp
    Sigma_im <- Sigma_i[start:end,start:end]
    yn[m] <- mu_i[m,jp+1] + Sigma_im[jp+1,index] %*% chol2inv(chol(Sigma_im[index,index])) %*% (yp[m,]-mu_i[m,index])
    mu_i_all[m,] <- theta_all[,m,] %*% X_tilde_i[jp+1,] + delta_all[,m,] %*% U_i[jp+1,] + alpha_i_all[,m,] %*% V_i[jp+1,]
  }
  
  # update the data 
  data$yp <- cbind(yp, yn)
  data$Xp <- X_i
  data$Hp <- H_i
  data$XHp <- XH_i
  data$Vp <- V_i
  data$tp <- t_i
  data$Zp <- Z_i
  data$Up <- U_i
  data$Dp <- D_i
  
  # uncertainty quantification
  uq <- apply(mu_i_all, 1, sd)
  
  # return list
  returnlist <- list(yn=yn, data=data, uq=uq)
  return(returnlist)
}



generate_state_pred <- function(data, jp, Zn, tnt, phi_nit){
  
  # generate the new state yn for prediction (uncertainty quantification not recorded) 
  # args: data: list of state, covarites, drug usuage, time at 1:jp-th visits
  #       Zn, tnt: regimen and time at jp+1-th visit
  #       phi_nit: posterior samples of state model parameters at one iteration
  # returns: yn: state matrix at jp+1-th visit
  #          data: the updated data matrix
  
  yn <- rep(NA, M)
  
  # data
  yp <- data$yp
  Xp <- data$Xp
  Hp <- data$Hp
  XHp <- data$XHp
  Vp <- data$Vp
  tp <- data$tp
  Zp <- data$Zp
  Up <- data$Up
  Dp <- data$Dp
  
  # parameters
  theta <- phi_nit$theta
  delta <- phi_nit$delta
  alpha_i <- phi_nit$alpha_i
  Cm <- phi_nit$Cm
  rho <- phi_nit$rho
  sigma2 <- phi_nit$sigma2
  
  # time correlation matrix
  t_i <- c(tp, tnt) # time points for i-th subject
  Ct_i <- auto_regressive_rcpp(n=jp+1, sigma2=1, rho, t_i)
  
  # drug toxicity
  Z_i <- rbind(Zp, c(Zn, rep(NA, max_num_drugs-length(Zn))))
  D_i <- rbind(Dp, as.numeric(drug_names %in% Zn))
  tit <- t_i - t_i[1] + 0.5; tmax <- max(tit)
  drug_toxic <- as.vector(drug_toxic_integrate_rcpp(num=100, lower=0, upper=tmax, Nk_sum, jp+1, D_i, tit))
  U_i <- rbind(Up, (drug_toxic-U_mean)/U_std)
  
  # fixed and random effect covariates
  X_i <- rbind(Xp, Xp[1,1:S]) 
  V_i <- rbind(Vp, c(1, (tnt-time_mean)/time_std))
  drug_simi <- Drug_Similarity(Zn, z[index_kernel,]) # drug similarity between Zn and kernel knots
  if (pca){
    # if pca on kernel regression 
    if (sum(drug_simi) != 0){ 
      drug_simi <- drug_simi/sum(drug_simi)
      drug_simi <- (drug_simi-H_mean)/H_std
      drug_simi <- drug_simi %*% eigenvec[,1:D] 
      H_i <- rbind(Hp, (drug_simi-Hpca_mean)/Hpca_std)
    }else{ H_i <- rbind(Hp, drug_simi[1:D]) }
  }else{
    # if no pca on kernel regression
    if (sum(drug_simi) != 0){ 
      drug_simi <- drug_simi/sum(drug_simi)
      H_i <- rbind(Hp, (drug_simi-H_mean)/H_std)
    }else{ H_i <- rbind(Hp, drug_simi) }
  }
  
  # drug comb and baseline interaction terms
  drug_base_inter <- Xp[1,2:S] %x% as.vector(drug_simi)
  XH_i <- rbind(XHp, (drug_base_inter-XH_mean)/XH_std)
  X_tilde_i <- cbind(X_i, H_i, XH_i)
  
  # mean and covariance matrix for Gaussian process
  mu_i <- matrix(NA, nrow=M, ncol=jp+1)
  Sigma_i <- Cm[1:M,1:M] %x% Ct_i + diag(sigma2, M) %x% diag(1, jp+1) 
  for (m in 1:M){
    # generate the new state matrix
    mu_i[m,] <- X_tilde_i %*% theta[m,] + U_i %*% delta[m,] + V_i %*% alpha_i[m,]
    start <- (m-1)*(jp+1)+1; end <- m*(jp+1); index <- 1:jp
    Sigma_im <- Sigma_i[start:end,start:end]
    yn[m] <- mu_i[m,jp+1] + Sigma_im[jp+1,index] %*% chol2inv(chol(Sigma_im[index,index])) %*% (yp[m,]-mu_i[m,index])
  }
  
  # update the data 
  data$yp <- cbind(yp, yn)
  data$Xp <- X_i
  data$Hp <- H_i
  data$XHp <- XH_i
  data$Vp <- V_i
  data$tp <- t_i
  data$Zp <- Z_i
  data$Up <- U_i
  data$Dp <- D_i

  # return list
  returnlist <- list(yn=yn, data=data)
  return(returnlist)
}



generate_data <- function(data, phi, theta, num_samples, num_pred){
  
  # generate predictive state and regimen data sequentially 
  # args: data: list of state, covarites, drug usuage, time at 1:Ji-th visits
  #       phi: posterior samples of state model parameters 
  #       theta: current regimen model parameters
  #       num_samples: number of posterior samples
  #       num_pred: number of forward predictions
  # returns: y_pred, Z_pred, uq_pred: predictive state matrices, regimens, and uncertainties
  
  y_pred <- list()
  Z_pred <- list()
  theta_all <- phi$theta[1:num_samples,,]
  delta_all <- phi$delta[1:num_samples,,]
  alpha_i_all <- phi$alpha[1:num_samples,i,,]
  uq_pred <- list()
  
  # generate data for each posterior samples
  for (nit in 1:num_samples){
    
    data_pred <- data
    data_pred$yp <- phi$y[nit,i,1:M,1:Ji]
    y_pred[[nit]] <- list()
    Z_pred[[nit]] <- list()
    uq_pred[[nit]] <- list()
    
    # state model nit-th posterior samples
    phi_nit <- NULL
    phi_nit$theta <- phi$theta[nit,,]
    phi_nit$delta <- phi$delta[nit,,]
    phi_nit$alpha_i <- phi$alpha[nit,i,,]
    phi_nit$Cm <- phi$Cm[nit,,]
    phi_nit$rho <- phi$rho[nit]
    phi_nit$sigma2 <- phi$sigma2[nit,]
    
    # generate a sequence of regimens and states 
    for (j in 1:num_pred){
      
      if (j == 1){
        # generate Z_{Ji+1} conditional on the observed state data y_{Ji}
        if (data_pred$yp[dep_index,Ji] <= dep_thres & 
            data_pred$yp[vload_index,Ji] <= vload_thres & 
            data_pred$yp[renal_index,Ji] >= renal_thres){
          # if not switch regimen 
          Z_pred[[nit]][[j]] <- data_pred$Zp[Ji,][!is.na(data_pred$Zp[Ji,])]
        }else{
          # if switch regimen 
          Z_pred[[nit]][[j]] <- generate_regimen(data_pred$yp[,Ji], theta)
        }
        # generate y_{Ji+1} conditional on the current regimen Z_{Ji+1}
        state_list <- generate_state(data_pred, Ji+j-1, Z_pred[[nit]][[j]], tn[j], phi_nit, theta_all, delta_all, alpha_i_all)
        y_pred[[nit]][[j]] <- state_list$yn; data_pred <- state_list$data; uq_pred[[nit]][[j]] <- state_list$uq
      }else{
        # generate Z_{Ji+j} conditional on the predictive state data y_{Ji+j-1}
        if (y_pred[[nit]][[j-1]][dep_index] <= dep_thres & 
            y_pred[[nit]][[j-1]][vload_index] <= vload_thres & 
            y_pred[[nit]][[j-1]][renal_index] >= renal_thres){
          # if not switch regimen 
          Z_pred[[nit]][[j]] <- Z_pred[[nit]][[j-1]]
        }else{
          # if switch regimen 
          Z_pred[[nit]][[j]] <- generate_regimen(y_pred[[nit]][[j-1]], theta)
        }
        # generate y_{Ji+1} conditional on the current regimen Z_{Ji+1}
        state_list <- generate_state(data_pred, Ji+j-1, Z_pred[[nit]][[j]], tn[j], phi_nit, theta_all, delta_all, alpha_i_all)
        y_pred[[nit]][[j]] <- state_list$yn; data_pred <- state_list$data; uq_pred[[nit]][[j]] <- state_list$uq
      }
    }
  }
  
  # returnlist 
  returnlist <- list(y_pred=y_pred, Z_pred=Z_pred, uq_pred=uq_pred)
  return(returnlist)
}



expected_reward <- function(y, Z, num_samples, num_pred, uq){
  
  # calculate the expected reward 
  # args: y: posterior predictive samples of state matrix
  #       Z: predictive ART regimens
  #       num_samples: number of posterior samples
  #       num_pred: number of forward predictions
  #       uq: uncertainty quantification
  # returns: r: rewards, R: expected reward, 
  
  r <- rep(NA, nrow=num_samples) # rewards
  for (nit in 1:num_samples){
    dep_reward <- (unlist(y[[nit]])[dep_index_all] - y_mean[dep_index]) / y_std[dep_index]
    vload_reward <- ((unlist(y[[nit]])[vload_index_all] > vload_thres) > 0) * 
      (abs(unlist(y[[nit]])[vload_index_all] - vload_thres) / y_std[vload_index])
    renal_reward <- ((unlist(y[[nit]])[renal_index_all] < renal_thres) > 0) * 
      (abs(unlist(y[[nit]])[renal_index_all] - renal_thres) / y_std[renal_index])
    r[nit] <- - sum(reward_coef * c(sum(dep_reward), sum(vload_reward), sum(renal_reward))) - 
      uq_factor*(sum(unlist(uq[[nit]])))
  }
  
  R <- mean(r) # expected_reward
  returnlist <- list(r=r, R=R)
  return(returnlist)
}



grad_expected_reward <- function(num_samples, num_pred, data, y, Z, r, R, theta){
  
  # calculate the gradient of expected reward 
  # args: num-samples: number of posterior samples
  #       num_pred: number of forward predictions
  #       data: list of state, covarites, drug usuage, time at 1:Ji-th visits
  #       y: posterior predictive samples of state matrix
  #       Z: posterior predictive samples of regimen
  #       theta: regimen model parameters
  # returns: grad: gradient of expected reward
  
  grad <- NULL
  grad$theta1 <- NULL
  grad$theta2 <- NULL
  
  # initialize gradients as zeros
  for (k in 1:K){
    grad$theta1[[k]] <- list()
    for (ck in 1:Ck[k]){ grad$theta1[[k]][[ck]] <- rep(0, M_tilde) }
    grad$theta2[[k]] <- list()
    for (nk in 1:Nk[k]){ grad$theta2[[k]][[nk]] <- rep(0, M_tilde) }
  }
  
  # calculate gradients 
  for (nit in 1:num_samples){
    data$yp <- phi$y[nit,i,1:M,1:Ji]
    for (j in 1:num_pred){
      # calculate gradient of log-likelihood of regimen model 
      if (j == 1){ grad_logll <- grad_logll_regimen(Z[[nit]][[j]], data$yp[,Ji], theta)
      }else{ grad_logll <- grad_logll_regimen(Z[[nit]][[j]], y[[nit]][[j-1]], theta) }
      re <- r[nit] - R # baseline subtraction
      # calculate gradient of reward function
      for (k in 1:K){
        for (ck in 1:Ck[k]){ grad$theta1[[k]][[ck]] <- grad$theta1[[k]][[ck]] + re*grad_logll$theta1[[k]][[ck]] }
        for (nk in 1:Nk[k]){ grad$theta2[[k]][[nk]] <- grad$theta2[[k]][[nk]] + re*grad_logll$theta2[[k]][[nk]] }
      }
    }
  }
  
  # take expectation on gradients
  for (k in 1:K){
    for (ck in 1:Ck[k]){ grad$theta1[[k]][[ck]] <- grad$theta1[[k]][[ck]]/num_samples }
    for (nk in 1:Nk[k]){ grad$theta2[[k]][[nk]] <- grad$theta2[[k]][[nk]]/num_samples }
  }
  
  return(grad)
}



initialize <- function(){
  
  # initialize SGD algorithm
  # returns: theta: regimen model parameters 
  #          y, Z: predictive state matrices and regimens
  #          r: rewards, R: expected reward
  
  # initialize theta
  theta <- list()
  for (k in 1:K){
    theta$theta1[[k]] <- list()
    for (ck in 1:Ck[k]){ theta$theta1[[k]][[ck]] <- rep(0, M_tilde) }
    theta$theta2[[k]] <- list()
    for (nk in 1:Nk[k]){ theta$theta2[[k]][[nk]] <- rep(0, M_tilde) }
  }
  
  # initialize states, regimens, and reward 
  data_list <- generate_data(data, phi, theta, num_samples, num_pred)
  y_pred <- data_list$y_pred; Z_pred <- data_list$Z_pred; uq_pred <- data_list$uq_pred
  reward_list <- expected_reward(y_pred, Z_pred, num_samples, num_pred, uq_pred)
  r <- reward_list$r; R <- reward_list$R
  
  # return init_list
  init_list <- list(theta=theta, y=y_pred, Z=Z_pred, r=r, R=R, uq=uq_pred)
  return(init_list)
}



update_theta <- function(theta, grad_theta, step_size){
  
  # update regimen model parameters theta using SGD
  # args: theta: regimen model parameters 
  #       grad_theta: gradients of theta
  #       step_size: step size of SGD
  
  theta_update <- theta
  
  # update theta using gradient descent
  for (k in 1:K){
    for (ck in 1:Ck[k]){ theta_update$theta1[[k]][[ck]] <- theta_update$theta1[[k]][[ck]] + step_size * grad_theta$theta1[[k]][[ck]] }
    for (nk in 1:Nk[k]){ theta_update$theta2[[k]][[nk]] <- theta_update$theta2[[k]][[nk]] + step_size * grad_theta$theta2[[k]][[nk]] }
  }
  
  return(theta_update)
}