####################################################### Data ########################################################

# drug class data 
K <- 6 # number of drug classes (two for NRTIs, booster RTV and COBI not count)
Ck <- rep(1, K) # number of maximum single drugs for each drug class in an ART regimen
max_num_drugs <- dim(Z)[3] # maximal number of drugs, i.e., dimension of the last column in array Z
M_tilde <- M + 1 # dimension for (add intercept) of y in the regimen model
all_drugs <- NULL # all possible single drugs for assignment in each drug class
all_drugs$NRTI1 <- c("3TC", "FTC") 
all_drugs$NRTI2 <- c("ABC", "TAF", "TDF") 
all_drugs$NNRTI <- c("EFV", "ETR", "NVP", "RPV")
all_drugs$PI <- c("ATV", "DRV", "LPV")
all_drugs$INSTI <- c("DTG", "EVG", "RAL")
all_drugs$EI <- c("MVC")
Nk <- rep(NA, K)  # number of all possible single drugs for assignment in each drug class
Nk <- c(length(all_drugs$NRTI1), length(all_drugs$NRTI2), length(all_drugs$NNRTI), 
        length(all_drugs$PI), length(all_drugs$INSTI), length(all_drugs$EI))

# reward function data 
dep_index <- 1; dep_index_all <- seq(dep_index, num_pred*M, by=M) # response index of depression
vload_index <- 2; vload_index_all <- seq(vload_index, num_pred*M, by=M) # response index of viral load
renal_index <- 3; renal_index_all <- seq(renal_index, num_pred*M, by=M) # response index of eGFR
dep_thres <- 16; vload_thres <- log(20); renal_thres <- 60 # thresholds for depression, viral load, and eGFR 
reward_coef <- reward_coef/sum(reward_coef) # reward function coefficients

# subject data 
Ji <- J[i] # number of visits 
pca <- TRUE # use PCA on drug combination kernel matrix
t_int <- 0.5 # time interval for perdictions
tn <- t[i,Ji] + seq(t_int, t_int*num_pred, by=t_int) # time points 
data <- NULL # subject data 
data$yp <- y[i,1:M,1:Ji]
data$Xp <- X[i,1:Ji,]; data$Hp <- H[i,1:Ji,]; data$XHp <- XH[i,1:Ji,] 
data$Vp <- V[i,1:Ji,]; data$tp <- t[i,1:Ji]; data$Zp <- Z[i,1:Ji,]
data$Up <- U[i,1:Ji,]; data$Dp <- drug_history[i,1:Ji,]
