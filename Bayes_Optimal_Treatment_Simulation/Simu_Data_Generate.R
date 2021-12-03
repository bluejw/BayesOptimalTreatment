
###################################### Load the Preprocessed Data ##################################

n <- data_preprocess$n # number of individuals
J <- data_preprocess$J # number of visits for each individual
Z <- data_preprocess$Z # ART regimens for each individual at each visit
z <- data_preprocess$z # unique ART regimens
kappa <- data_preprocess$kappa # subset-tree kernel similarity matrix
index_regimens <- data_preprocess$index_regimens # index of ART regimens in the unique regimen array
all_drugs <- data_preprocess$all_drugs # list of all individual ART drugs
abbr_group <- data_preprocess$abbr_group # drug class data
drug_names <- as.vector(unlist(all_drugs)) # single ART drug names sorted by drug class
K <- 6 # number of drug classes, including NRTI, NNRTI, PI, INSTI, EI, PE
Nk <- rep(NA, K) # number of all possible single drugs for each drug class
Nk <- c(length(all_drugs$NRTI), length(all_drugs$NNRTI), length(all_drugs$PI), 
        length(all_drugs$INSTI), length(all_drugs$EI), length(all_drugs$PE)) 
Nk_sum <- sum(Nk) # number of all single ART drugs

##################################### Start to Generate Simulated Data ###############################

set.seed(123) # set seed 

### basline covariates 
S <- 3 # dimension of basline covariates 
X <- array(NA, dim=c(n, max(J), S)) # the covaraite matrix
for (i in 1:n){
  X[i,1:J[i],1] <- rep(1, J[i]) # intercept term
  X[i,1:J[i],2] <- rep(sample(0:1, 1, prob=c(0.6,0.4), replace=T), J[i]) # time-invariant binary covariate
  X[i,1:J[i],3] <- rep(rnorm(1, 0, 1), J[i]) # time-invariant continuous covariate
}

# normalize covariates X for each dimension
X_mean <- rep(NA, S); X_std <- rep(NA, S)
for (s in 2:S){
  X_hat <- X[,,s][!is.na(X[,,s])]
  X[,,s] <- (X[,,s] - mean(X_hat)) / sd(X_hat)
  X_mean[s] <- mean(X_hat); X_std[s] <- sd(X_hat)
}

### drug combinations
drug_times <- as.vector(table(index_regimens)) # times of drug apperance in data
drug_index_all <- sort(unique(as.vector(index_regimens))) # index of all drug
index_kernel <- drug_index_all[which(drug_times>=10)] # index of representative ART regimens
D <- length(index_kernel) # number of drug used in kernel regression: D = 67

H <- array(NA, dim=c(n, max(J), D)) # drug kernel regression design matrix
for (i in 1:n){ 
  for (j in 1:J[i]){
    if (sum(kappa[index_regimens[i,j],index_kernel])!=0){
      H[i,j,1:D] <- kappa[index_regimens[i,j],index_kernel] / sum(kappa[index_regimens[i,j],index_kernel]) 
    }else{
      H[i,j,1:D] <- 0 
    }
  }
}

# normalize covariates H for each dimension
H_mean <- rep(NA, D); H_std <- rep(NA, D)
for (d in 1:D){
  H_hat <- H[,,d][!is.na(H[,,d])]
  H[,,d] <- (H[,,d] - mean(H_hat)) / sd(H_hat)
  H_mean[d] <- mean(H_hat); H_std[d] <- sd(H_hat)
}

### PCA
index_start <- 1; index_end <- sum(J)
H_tilde <- matrix(NA, nrow=index_end, ncol=D) # sum(J)*D matrix 
for (i in 1:n){
  H_tilde[index_start:(index_start+J[i]-1),] <- H[i,1:J[i],]
  index_start <- index_start + J[i]
}
  
H_svd <- svd(H_tilde)
eigenval <- H_svd$d^2
sum(eigenval[1:41])/sum(eigenval); plot(eigenval) # retain 99.9% variation 
eigenvec <- H_svd$v; D_pca <- 41
H_tilde_proj <- H_tilde %*% eigenvec[,1:D_pca] 
  
index_start <- 1; index_end <- sum(J)
H_proj <- array(NA, dim=c(n, max(J), D_pca))
for (i in 1:n){
  H_proj[i,1:J[i],] <- H_tilde_proj[index_start:(index_start+J[i]-1),]
  index_start <- index_start + J[i]
}
  
D <- D_pca; H <- H_proj
Hpca_mean <- rep(NA, D); Hpca_std <- rep(NA, D)
for (d in 1:D){
  H_hat <- H[,,d][!is.na(H[,,d])]
  H[,,d] <- (H[,,d] - mean(H_hat)) / sd(H_hat)
  Hpca_mean[d] <- mean(H_hat); Hpca_std[d] <- sd(H_hat)
}

### baseline and drug combination interation terms
XH <- array(NA, dim=c(n, max(J), (S-1)*D))
for (i in 1:n){
  for (j in 1:J[i]){
    XH[i,j,] <- X[i,j,2:S] %x% H[i,j,]
  }
}

# normalize covariates XH for each dimension
XH_mean <- rep(NA, (S-1)*D); XH_std <- rep(NA, (S-1)*D)
for (sd in 1:((S-1)*D)){
  XH_hat <- XH[,,sd][!is.na(XH[,,sd])]
  XH[,,sd] <- (XH[,,sd] - mean(XH_hat)) / sd(XH_hat)
  XH_mean[sd] <- mean(XH_hat); XH_std[sd] <- sd(XH_hat)
}

### random effect covariates 
Q <- 2 # dimension of random effect covariates 
V <- array(NA, dim=c(n, max(J), Q)) # random effect covariates
time_mean <- 0; time_std <- 1 # mean and std of time
t <- matrix(NA, nrow=n, ncol=max(J)) # time 
for (i in 1:n){
  t[i,1:J[i]] <- sort(rnorm(J[i], time_mean, time_std))
  V[i,1:J[i],1] <- rep(1, J[i]) # intercept term
  V[i,1:J[i],2] <- t[i,1:J[i]] # time 
}

### drug toxicity 
drug_history <- array(NA, dim=c(n, max(J), Nk_sum)) # drug history indicator matrix
for (i in 1:n){ for (j in 1:J[i]){ drug_history[i,j,] <- as.numeric(drug_names %in% Z[i,j,]) }}
ti <- matrix(NA, nrow=n, ncol=max(J)) # time for each subject used for calculating drug toxicity
for (i in 1:n){ ti[i,1:J[i]] <- t[i,1:J[i]] - t[i,1] + 0.5 }
U <- array(NA, dim=c(n, max(J), Nk_sum)) # drug toxicity covariates
for (i in 1:n){
  for (j in 1:J[i]){
    if (j == 1){ drug_history_ij <- as.matrix(t(drug_history[i,1,]))
    }else{ drug_history_ij <- drug_history[i,1:j,] }
    U[i,j,] <- drug_toxic_integrate_rcpp(num=100, lower=0, upper=ti[i,j], Nk_sum, j, drug_history_ij, ti[i,1:j])
  }
}

# normalize covariates U for each dimension
U_mean <- rep(NA, Nk_sum); U_std <- rep(NA, Nk_sum)
for (nk in 1:Nk_sum){
  U_hat <- U[,,nk][!is.na(U[,,nk])]
  U[,,nk] <- (U[,,nk] - mean(U_hat)) / sd(U_hat)
  U_mean[nk] <- mean(U_hat); U_std[nk] <- sd(U_hat)
}

##################################### Simulated True Parameters ###################################

M <- 3 # dimension of states: depression, viral load, renal function

### fixed effect coefficients
beta <- matrix(NA, nrow=M, ncol=S) # baseline coefficients
beta[1,] <- c(25, 1, 2); beta[2,] <- c(4.5, -0.5, 1); beta[3,] <- c(75, -4, 2)
gamma <- matrix(NA, nrow=M, ncol=D) # drug combination coefficients
for (m in 1:M){ gamma[m,] <- rmvn_rcpp(1, rep(0,D), diag(1,D)) }
phi <- matrix(0, nrow=M, ncol=(S-1)*D) # baseline and drug combination interation coefficients
for (m in 1:M){ phi[m,] <- rmvn_rcpp(1, rep(0,(S-1)*D), diag(1,(S-1)*D)) }
delta <- matrix(0, nrow=M, ncol=Nk_sum) # drug toxicity coefficients
delta[1,1:Nk[1]] <- rep(1, Nk[1]); delta[2,1:Nk[1]] <- rep(0.5, Nk[1]); delta[3,1:Nk[1]] <- rep(-2, Nk[1])

### random effect coefficients
Sigma_alpha <- array(NA, dim=c(M, Q, Q)) # random effect matrix 
Sigma_alpha[1,,] <- diag(1, Q); Sigma_alpha[2,,] <- diag(0.5, Q); Sigma_alpha[3,,] <- diag(2, Q)
alpha <- array(NA, dim=c(n, M, Q)) # random effect coefficients
for (i in 1:n){ for (m in 1:M){ alpha[i,m,] <- rmvn_rcpp(1, rep(0,Q), Sigma_alpha[m,,]) }}
alpha[1,2,] <- c(-4, -1); alpha[2,2,] <- c(-2.5, -1); alpha[2,3,] <- c(-8, -2) # subject #1 and #2 

### GP parameters
Cm <- matrix(1, nrow=M, ncol=M) # state GP correlation matrix 
Cm[1,2] <- Cm[2,1] <- 0.75; Cm[1,3] <- Cm[3,1] <- -0.5; Cm[2,3] <- Cm[3,2] <- -0.25
sigma2_omega <- c(5, 1, 10) # variance of GP correlation terms
Cm <- (diag(sqrt(sigma2_omega))%*%Cm%*%diag(sqrt(sigma2_omega)))
rho <- 0.5 # GP time correlation
Ct <- array(NA, dim=c(n, max(J), max(J))) # GP time correlation matrix
for (i in 1:n){ for (j in 1:J[i]){ Ct[i,j,1:J[i]] <- rho^(abs(t[i,j]-t[i,1:J[i]])) }}
omega <- array(NA, dim=c(n, M, max(J))) # GP correlation terms
for (i in 1:n){
  omega_i <- rmvn_rcpp(1, rep(0,M*J[i]), Cm %x% Ct[i,1:J[i],1:J[i]])
  omega[i,1:M,1:J[i]] <- matrix(omega_i, nrow=M, ncol=J[i], byrow=T)
}

### simulated true state matrix 
sigma2 <- c(10, 1, 20) # i.i.d error variance 
y <- array(NA, dim=c(n, M, max(J))) # state matrix 
data_index <- array(0, dim=c(n, M, max(J))) # index of non-missing data for state matrix
for (i in 1:n){
  for (m in 1:M){
    mu_im <- X[i,1:J[i],] %*% beta[m,] + H[i,1:J[i],] %*% gamma[m,] + XH[i,1:J[i],] %*% phi[m,] +
      U[i,1:J[i],] %*% delta[m,] + V[i,1:J[i],] %*% alpha[i,m,] + omega[i,m,1:J[i]]
    # generate index of non-missing data 
    num_obs <- sample(1:J[i], 1, replace = FALSE, prob=dpois(1:J[i], lambda=25))
    data_index_im <- sort(sample(1:J[i], num_obs, replace = FALSE))
    data_index[i,m,data_index_im] <- rep(1, num_obs)
    # generate simulated true state matrix
    y[i,m,data_index_im] <- rmvn_rcpp(1, mu_im[data_index_im], diag(sigma2[m],num_obs))
  }
}

y_mean <- rep(NA, M); y_std <- rep(NA, M)
y_tilde <- array(NA, dim=c(n, M, max(J)))
for (m in 1:M){
  y_hat <- y[,m,][!is.na(y[,m,])]
  y_tilde[,m,] <- (y[,m,] - mean(y_hat)) / sd(y_hat)
  y_mean[m] <- mean(y_hat); y_std[m] <- sd(y_hat)
}
y_min <- min(y_tilde[!is.na(y_tilde)])
y_max <- max(y_tilde[!is.na(y_tilde)])

###################################### Hyper-parameters and Pre-calculated data ##################################

### prior specification: hyper-parameters
S_tilde <- S+D+(S-1)*D # dimension of fixed effects
shrink_index <- 2:S_tilde # index of shrinkage parameters
S_minus <- length(shrink_index) # dimension of shrinkage parameters
theta <- cbind(beta, gamma, phi) # S_tilde-dimensional fixed effects 
mu_theta <- rep(0, S_tilde); Sigma_theta <- diag(100, S_tilde) # fixed effects prior
Sigma_theta_inv <- chol2inv(chol(Sigma_theta)); Sigma_theta_inv_mu_theta <- Sigma_theta_inv %*% mu_theta
a_0 <- Q+1; A_0 <- diag(1, Q) # Sigma_alpha \sim Inverse-Wishart(a_0, A_0^{-1}) # random effects prior 
g_1 <- 1; g_2 <- 1 # sigma2_m \sim inverse-gamma(g_1, g_2)
L <- 1 # dimension of approximation matrix for GP covariance matrix
negative_state <- c(3) # negative state, e.g., CD4, EGFR

### pre-calculated data for saving running time  
X_tilde <- array(NA, dim=c(n, max(J), S_tilde))
for (i in 1:n){ X_tilde[i,1:J[i],] <- cbind(X[i,1:J[i],], H[i,1:J[i],], XH[i,1:J[i],]) }

######################################## Save Simulation Truths ######################################

data <- NULL # simulation truths

# data
data$S <- S
data$X <- X; data$X_mean <- X_mean; data$X_std <- X_std
data$index_kernel <- index_kernel
data$D <- D
data$eigenvec <- eigenvec
data$H <- H; data$H_mean <- H_mean; data$H_std <- H_std
data$XH  <- XH; data$XH_mean <- XH_mean; data$XH_std <- XH_std
data$Q <- Q
data$V <- V
data$t <- t; data$time_mean <- time_mean; data$time_std <- time_std
data$drug_history <- drug_history 
data$U <- U; data$U_mean <- U_mean; data$U_std <- U_std
data$M <- M
data$data_index <- data_index
data$y <- y; data$y_mean <- y_mean; data$y_std <- y_std

# parameters
data$beta <- beta
data$gamma <- gamma
data$phi <- phi
data$delta <- delta
data$alpha <- alpha
data$Sigma_alpha <- Sigma_alpha
data$Cm <- Cm
data$rho <- rho
data$omega <- omega
data$sigma2 <- sigma2

save(data, file="Simu_Truths.Rdata")