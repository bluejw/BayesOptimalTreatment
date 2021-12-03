####################################################### SGD ########################################################


SGD_Main <- function(i, phi, Niter, num_samples, num_pred, step_size, reward_coef, uq_factor){
  
  # SGD main function 
  # args: i: index of subject, phi: posterior samples 
  #       Niter: number of SGD iterations
  #       num_samples: number of posterior samples
  #       num_pred: number of forward predictions
  #       step_size: step size of SGD
  #       reward_coef: reward weight for subject i
  #       uq_factor: uncertainty quantification factor
  # returns: SGD results 
  
  ### SGD
  sgd <- NULL
  sgd$theta <- NULL
  sgd$grad <- NULL
  sgd$y <- NULL
  sgd$Z <- NULL
  sgd$r <- NULL
  sgd$R <- NULL
  
  ### Initialize
  init <- initialize()
  sgd$theta[[1]] <- init$theta
  sgd$grad[[1]] <- list()
  sgd$y[[1]] <- init$y
  sgd$Z[[1]] <- init$Z
  sgd$r[[1]] <- init$r
  sgd$R[[1]] <- init$R
  sgd$uq[[1]] <- init$uq
  
  ### Start of the algorithm
  
  start.time = proc.time()
  
  for (iter in 2:Niter){
    
    print(iter)
    
    # generate state matrix and regimen 
    data_list <- generate_data(data, phi, sgd$theta[[iter-1]], num_samples, num_pred)
    sgd$y[[iter]] <- data_list$y_pred; sgd$Z[[iter]] <- data_list$Z_pred; sgd$uq[[iter]] <- data_list$uq_pred
    reward_list <- expected_reward(sgd$y[[iter]], sgd$Z[[iter]], num_samples, num_pred, sgd$uq[[iter]])
    sgd$r[[iter]] <- reward_list$r; sgd$R[[iter]] <- reward_list$R
    
    # calculate the gradient of expected reward
    sgd$grad[[iter]] <- grad_expected_reward(num_samples, num_pred, data, 
                sgd$y[[iter]], sgd$Z[[iter]], sgd$r[[iter]], sgd$R[[iter]], sgd$theta[[iter-1]])
    
    # update the parameters by gradient descent
    sgd$theta[[iter]] <- update_theta(sgd$theta[[iter-1]], sgd$grad[[iter]], step_size)
  }  
  
  duration = proc.time()-start.time
  print(duration) # print the running time 
  
  ### End of the algorithm
  
  return(sgd)
}
