
###################################### Load the Preprocessed Data ##################################

n <- data_preprocess$n # number of individuals
J <- data_preprocess$J # number of visits for each individual
Z <- data_preprocess$Z # ART regimens for each individual at each visit
z <- data_preprocess$z # unique ART regimens
kappa <- data_preprocess$kappa # subset-tree kernel similarity matrix
index_regimens <- data_preprocess$index_regimens # index of ART regimens in the unique regimen array
abbr_group <- data_preprocess$abbr_group # drug class data
combined_data <- data_preprocess$combined_data # combined data 

######################################### State Matrix ############################################

y_index <- c("DEPRESSION", "VLOAD", "EGFR2", "BMI") # index of states
M <- length(y_index) # number of states 
y <- array(NA, dim=c(n, M, max(J))) # state matrix 
data_index <- array(0, dim=c(n, M, max(J))) # index of non-missing data for state matrix
for (i in 1:n){
  # data for i-th individual
  current_data <- combined_data[which(combined_data$ID==i),] 
  data_index[i,1:M,1:J[i]] <- t(!is.na(current_data[,y_index]))
  y[i,1,1:J[i]] <- current_data$DEPRESSION
  y[i,2,1:J[i]] <- log(current_data$VLOAD) 
  y[i,3,1:J[i]] <- current_data$EGFR2
  y[i,4,1:J[i]] <- current_data$BMI
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

################################## Fixed Effects: Baseline #######################################

x_index <- c("INTERCEPT", "AGEATVIS", "DIABETES", "CURSMOKE", "DRUGUSE", "ANYHTN", "EMPLOY")
S <- length(x_index) # dimension of baseline covariates 
X <- array(NA, dim=c(n, max(J), S)) # the covaraite matrix
for (i in 1:n){
  current_data <- combined_data[which(combined_data$ID==i),]
  X[i,1:J[i],1] <- rep(1, J[i])
  X[i,1:J[i],2] <- rep(current_data$AGEATVIS[1], J[i])
  X[i,1:J[i],3] <- rep(current_data$DIABETES[1], J[i])
  X[i,1:J[i],4] <- rep(current_data$CURSMOKE[1], J[i])
  X[i,1:J[i],5] <- rep(current_data$DRUGUSE[1], J[i])
  X[i,1:J[i],6] <- rep(current_data$ANYHTN[1], J[i])
  X[i,1:J[i],7] <- rep(current_data$EMPLOY[1], J[i])
}

# normalize covariates X for each dimension
for (s in c(2:S)){
  X_hat <- X[,,s][!is.na(X[,,s])]
  X[,,s] <- (X[,,s] - mean(X_hat)) / sd(X_hat)
}

############################### Fixed Effects: Drug Combination ####################################

drug_times <- as.vector(table(index_regimens)) # times of drug apperance in data
drug_index_all <- sort(unique(as.vector(index_regimens))) # index of all drug
index_kernel <- drug_index_all[which(drug_times>=10)]
D <- length(index_kernel) # number of drug used in kernel regression: D = 105

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

# PCA
index_start <- 1; index_end <- sum(J)
H_tilde <- matrix(NA, nrow=index_end, ncol=D) # sum(J)*D matrix 
for (i in 1:n){
  H_tilde[index_start:(index_start+J[i]-1),] <- H[i,1:J[i],]
  index_start <- index_start + J[i]
}
  
H_svd <- svd(H_tilde)
eigenval <- H_svd$d^2
sum(eigenval[1:51])/sum(eigenval); plot(eigenval) # retain 99.9% variation 
eigenvec <- H_svd$v; D_pca <- 51
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

######################################## Random Effects ############################################

t <- matrix(NA, nrow=n, ncol=max(J)) # time 
first_date <- min(as.Date(combined_data$DATE)) # first date in the combined data
for (i in 1:n){
  time_i <- combined_data$DATE[which(combined_data$ID==i)]
  t[i,1:J[i]] <- as.integer(difftime(time_i, first_date, units = c("days"))) / 365
}

Q <- 2 # dimension of random effect covariates 
V <- array(NA, dim=c(n, max(J), Q)) # random effect covariates
for (i in 1:n){
  V[i,1:J[i],1] <- rep(1, J[i]) # intercept term
  V[i,1:J[i],2] <- t[i,1:J[i]] # time 
}

# normalize covariates V for each dimension
V_hat <- V[,,2][!is.na(V[,,2])]
V[,,2] <- (V[,,2] - mean(V_hat)) / sd(V_hat)
time_mean <- mean(V_hat); time_std <- sd(V_hat)

################################## Fixed Effects: Drug Toxicity ######################################

drug_names_data <- sort(unique(as.vector(Z[!is.na(Z)]))) # single ART drug names in data 
abbr_group <- abbr_group[which(names(abbr_group) %in% drug_names_data)]
abbr_group <- abbr_group[,order(names(abbr_group))] 

all_drugs <- NULL # all possible single drugs for each drug class
all_drugs$NRTI <- sort(drug_names_data[which(abbr_group=="NRTI",)])
all_drugs$NNRTI <- sort(drug_names_data[which(abbr_group=="NNRTI",)])
all_drugs$PI <- sort(drug_names_data[which(abbr_group=="PI",)])
all_drugs$INSTI <- sort(drug_names_data[which(abbr_group=="INSTI",)])
all_drugs$EI <- sort(drug_names_data[which(abbr_group=="EI",)])
all_drugs$PE <- sort(drug_names_data[which(abbr_group=="PE",)])
drug_names <- as.vector(unlist(all_drugs)) # single ART drug names sorted by drug class

K <- 6 # number of drug classes, including NRTI, NNRTI, PI, INSTI, EI, PE
Nk <- rep(NA, K) # number of all possible single drugs for each drug class
Nk <- c(length(all_drugs$NRTI), length(all_drugs$NNRTI), length(all_drugs$PI), 
        length(all_drugs$INSTI), length(all_drugs$EI), length(all_drugs$PE))
Nk_sum <- sum(Nk) # number of all single ART drugs

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

#################################### Fixed Effects: Interations ######################################

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

########################################### Computation Setup ##########################################

# prior specification: hyper-parameters
S_tilde <- S+D+(S-1)*D # dimension of fixed effects
shrink_index <- 2:S_tilde # index of shrinkage parameters
S_minus <- length(shrink_index) # dimension of shrinkage parameters
a_0 <- Q+1; A_0 <- diag(1, Q) # Sigma_alpha \sim Inverse-Wishart(a_0, A_0^{-1}) # random effects prior 
g_1 <- 1; g_2 <- 1 # sigma2_m \sim inverse-gamma(g_1, g_2)
L <- 1 # dimension of approximation matrix for GP covariance matrix
negative_state <- which(y_index %in% c("EGFR2", "BMI")) # negative state, e.g., eGFR, BMI

# pre-calculated data for saving running time  
X_tilde <- array(NA, dim=c(n, max(J), S_tilde))
for (i in 1:n){ X_tilde[i,1:J[i],] <- cbind(X[i,1:J[i],], H[i,1:J[i],], XH[i,1:J[i],]) }