##########################################################################################
#                                  WIHS Data Analysis                                    #
##########################################################################################

library(Rcpp)
library(RcppArmadillo)
library(MCMCpack)
library(coda)
library(TruncatedNormal)
library(Matrix)
library(BiasedUrn)
library(ggplot2)
library(gridExtra)
library(viridis)

source("MCMC_R_Functions.R") # R functions for MCMC 
sourceCpp("MCMC_Rcpp_Functions.cpp") # Rcpp functions for MCMC 
load("WIHS_Data_Preprocess.Rdata") # load the WIHS data 
source("WIHS_Data_Preprocess.R") # preprocess WIHS data 

##########################################################################################
#                           MCMC Posterior Inference                                     #
##########################################################################################

load("WIHS_MCMC_Results.Rdata") # posterior samples

mean <- NULL # posterior means of parameters 
mean$theta <- colMeans(post$theta)
mean$delta <- colMeans(post$delta)

Cred_Interval <- function(post_samples){
  # 95% credible intervals for posterior samples of matrix 
  num_row <- dim(post_samples)[2]; num_col <- dim(post_samples)[3]
  cred_interval <- array(NA, dim=c(num_row, num_col, 2)) # credible intervals for matrix 
  for (i in 1:num_row){ for (j in 1:num_col){ cred_interval[i,j,] <- quantile(post_samples[,i,j], c(0.025, 0.975)) }}
  return(cred_interval)
}

ci <- NULL # 95% credible intervals of parameters 
ci$theta <- Cred_Interval(post$theta)
ci$delta <- Cred_Interval(post$delta)

##########################################################################################
#         Figure 6: Posterior Means and 95% CIs for the Estimated Coefficients           #
##########################################################################################

# (a) Baseline coefficients
df <- data.frame(par = rep(c("Age","Smoke","Hypertension","Employment"),4), value = as.vector(t(mean$theta[,c(2,4,6,7)])),
                 lower = as.vector(t(ci$theta[,c(2,4,6,7),1])), upper = as.vector(t(ci$theta[,c(2,4,6,7),2])),
                 state = factor(c(rep("Depression",4), rep("Viral Load",4), rep("eGFR",4), rep("BMI",4)),
                                levels = c("Depression","Viral Load","eGFR","BMI")))
p <- ggplot(df, aes(par, value, colour = state)) + coord_flip() + 
  theme_classic() + geom_hline(yintercept=0, colour="darkgray") + 
  geom_pointrange(aes(ymin = lower, ymax = upper), position=position_dodge(0.5), size=1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0, position=position_dodge(0.5), size=1.5) + 
  labs(caption = " ", x = NULL, y = NULL, col = " ") + 
  theme(legend.title=element_text(size = 25, face='bold')) + theme(legend.text=element_text(size = 25)) + 
  theme(axis.title.x = element_text(size = 25, face='bold')) + theme(axis.title.y = element_text(size = 25, face='bold')) +
  theme(axis.text.x = element_text(size = 25)) + theme(axis.text.y = element_text(size = 25))

name <- paste("real_esti_base.pdf")
pdf(name, width = 18, height = 6, onefile = TRUE)
p
dev.off()

# (b) Drug combination coefficients
df <- data.frame(par = rep(c("PC1","PC2","PC5","PC6"),4), value = as.vector(t(mean$theta[,c(8,9,12,13)])),
                 lower = as.vector(t(ci$theta[,c(8,9,12,13),1])), upper = as.vector(t(ci$theta[,c(8,9,12,13),2])),
                 state = factor(c(rep("Depression",4), rep("Viral Load",4), rep("eGFR",4), rep("BMI",4)),
                                levels = c("Depression","Viral Load","eGFR","BMI")))
p <- ggplot(df, aes(par, value, colour = state)) + coord_flip() + 
  theme_classic() + geom_hline(yintercept=0, colour="darkgray") + 
  geom_pointrange(aes(ymin = lower, ymax = upper), position=position_dodge(0.5), size=1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0, position=position_dodge(0.5), size=1.5) + 
  labs(caption = " ", x = NULL, y = NULL, col = " ") + 
  theme(legend.title=element_text(size = 25, face='bold')) + theme(legend.text=element_text(size = 25)) + 
  theme(axis.title.x = element_text(size = 25, face='bold')) + theme(axis.title.y = element_text(size = 25, face='bold')) +
  theme(axis.text.x = element_text(size = 25)) + theme(axis.text.y = element_text(size = 25))

name <- paste("real_esti_drugcomb.pdf")
pdf(name, width = 18, height = 6, onefile = TRUE)
p
dev.off()

# (c) Drug toxicity coefficients
df <- data.frame(par = rep(c("D4T","TDF","DLV","MVC"),4), value = as.vector(t(mean$delta[,c(4,9,10,28)])),
                 lower = as.vector(t(ci$delta[,c(4,9,10,28),1])), upper = as.vector(t(ci$delta[,c(4,9,10,28),2])),
                 state = factor(c(rep("Depression",4), rep("Viral Load",4), rep("eGFR",4), rep("BMI",4)),
                                levels = c("Depression","Viral Load","eGFR","BMI")))
p <- ggplot(df, aes(par, value, colour = state)) + coord_flip() + 
  theme_classic() + geom_hline(yintercept=0, colour="darkgray") + 
  geom_pointrange(aes(ymin = lower, ymax = upper), position=position_dodge(0.5), size=1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0, position=position_dodge(0.5), size=1.5) + 
  labs(caption = " ", x = NULL, y = NULL, col = " ") + 
  theme(legend.title=element_text(size = 25, face='bold')) + theme(legend.text=element_text(size = 25)) + 
  theme(axis.title.x = element_text(size = 25, face='bold')) + theme(axis.title.y = element_text(size = 25, face='bold')) +
  theme(axis.text.x = element_text(size = 25)) + theme(axis.text.y = element_text(size = 25))

name <- paste("real_esti_drugtoxic.pdf")
pdf(name, width = 18, height = 6, onefile = TRUE)
p
dev.off()

##########################################################################################
#        Table 1: Top Three Positively and Negatively Related cART Regimens              #
##########################################################################################

# the index of top three positively and negatively related
index <- c(105,104,103,1,2,3) 

# PC1 
z[index_kernel[order(eigenvec[,1])[index]],]
sort(eigenvec[,1])[index]
# PC2
z[index_kernel[order(eigenvec[,2])[index]],]
sort(eigenvec[,2])[index]
# PC5
z[index_kernel[order(eigenvec[,5])[index]],]
sort(eigenvec[,5])[index]
# PC6
z[index_kernel[order(eigenvec[,6])[index]],]
sort(eigenvec[,6])[index]

##########################################################################################
#                   Figure 7: SGD Results (uncertainty is not penalized)                 #
##########################################################################################

source("SGD_R_Functions.R") # R functions for SGD
sourceCpp("SGD_Rcpp_Functions.cpp") # Rcpp functions for SGD
source("Drug_Similarity.R") # subset-tree drug similarity 
source("SGD_Main.R") # SGD main function
load("WIHS_SGD_Results.Rdata") # SGD results for WIHS data analysis

# Individual S1

i <- 3 # subject index 
phi <- post # posterior samples 
Niter <- 1000 # number of SGD iterations
num_samples <- 100 # number of posterior samples
num_pred <- 4 # number of forward predictions
step_size <- 0.1 # step size of gradient descent
reward_coef <- c(0.4,0.5,0.1) # reward weights
uq_factor <- 0 # uncertainty quantification factors
source("SGD_Data.R") # Data processing for SGD 
sgd <- sgd_all$S1_uq0 

# (a) Expected reward for S1
df <- data.frame(iteration=seq(1,1000,by=1), reward=unlist(sgd$R))
p <- ggplot(df, aes(iteration, reward)) + geom_line(colour = "black") + 
  theme_classic() + xlab("SGD Iteration") + ylab("Expected Reward") + 
  theme(axis.title.x = element_text(size = 30)) + theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30)) 

name <- paste("real_sgd_reward_1.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

# (c) Viral load for S1
y_pred <- list(); Z_pred <- list()
yt <- array(NA, dim=c(num_samples, num_pred, M))
# generate data for each posterior samples
for (nit in 1:num_samples){
  
  data_pred <- data
  data_pred$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred[[nit]] <- list(); Z_pred[[nit]] <- list()
  Z_pred[[nit]][[1]] <- c("ABC","FTC","RAL")
  Z_pred[[nit]][[2]] <- c("ABC","FTC","MVC","RAL")
  Z_pred[[nit]][[3]] <- c("ABC","FTC","MVC","RAL")
  Z_pred[[nit]][[4]] <- c("ABC","FTC","MVC","RAL")
  
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
    # generate y_{Ji+1} conditional on the current regimen Z_{Ji+1}
    state_list <- generate_state_pred(data_pred, Ji+j-1, Z_pred[[nit]][[j]], tn[j], phi_nit)
    y_pred[[nit]][[j]] <- state_list$yn; data_pred <- state_list$data; yt[nit,j,] <- y_pred[[nit]][[j]]
  }
}

y_vload <- data$yp[vload_index,]
yt_mean_vload <- colMeans(yt[,,vload_index])
yt_lower_vload <- apply(yt[,,vload_index], 2, function(x){quantile(x, 0.025)})
yt_upper_vload <- apply(yt[,,vload_index], 2, function(x){quantile(x, 0.975)})

df <- data.frame(x=seq(1,Ji+num_pred,by=1), y=c(y_vload, yt_mean_vload),
                 lower=c(rep(NA,Ji),yt_lower_vload), upper=c(rep(NA,Ji),yt_upper_vload))
p <- ggplot(df, aes(x, y)) + 
  geom_ribbon(aes(x, ymin=lower, ymax=upper), fill="gray", colour="gray") + 
  geom_line(size=2) + geom_hline(yintercept=log(20), colour='red', linetype='dashed') +
  theme_classic() + xlab("Visit") + ylab("Viral Load") + 
  theme(axis.title.x = element_text(size = 30)) + theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30)) +
  scale_x_continuous(labels = as.character(1:length(df$x)), breaks = 1:length(df$x))

name <- paste("real_state_pred_1.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

# Individual S2

i <- 84 # subject index 
reward_coef <- c(1/3,1/3,1/3) # reward weights
uq_factor <- 0 # uncertainty quantification factor
source("SGD_Data.R") # Data processing for SGD 
sgd <- sgd_all$S2_uq0

# (b) Expected reward for S2
df <- data.frame(iteration=seq(1,1000,by=1), reward=unlist(sgd$R))
p <- ggplot(df, aes(iteration, reward)) + geom_line(colour = "black") +  
  theme_classic() + xlab("SGD Iteration") + ylab("Expected Reward") + 
  theme(axis.title.x = element_text(size = 30)) + theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30)) 

name <- paste("real_sgd_reward_2.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

# (d) Depression for S2
y_pred <- list(); Z_pred <- list()
yt <- array(NA, dim=c(num_samples, num_pred, M))
# generate data for each posterior samples
for (nit in 1:num_samples){
  
  data_pred <- data
  data_pred$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred[[nit]] <- list(); Z_pred[[nit]] <- list()
  Z_pred[[nit]][[1]] <- c("ABC","ATV","DTG","FTC","RTV")
  Z_pred[[nit]][[2]] <- c("ABC","ATV","COBI","EVG","FTC","RTV")
  Z_pred[[nit]][[3]] <- c("ABC","ATV","COBI","EVG","FTC","RTV")
  Z_pred[[nit]][[4]] <- c("ABC","COBI","EVG","FTC","LPV","RTV")
  
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
    # generate y_{Ji+1} conditional on the current regimen Z_{Ji+1}
    state_list <- generate_state_pred(data_pred, Ji+j-1, Z_pred[[nit]][[j]], tn[j], phi_nit)
    y_pred[[nit]][[j]] <- state_list$yn; data_pred <- state_list$data; yt[nit,j,] <- y_pred[[nit]][[j]]
  }
}

y_dep <- data$yp[dep_index,]
yt_mean_dep <- colMeans(yt[,,dep_index])
yt_lower_dep <- apply(yt[,,dep_index], 2, function(x){quantile(x, 0.025)})
yt_upper_dep <- apply(yt[,,dep_index], 2, function(x){quantile(x, 0.975)})

df <- data.frame(x=seq(1,Ji+num_pred,by=1), y=c(y_dep, yt_mean_dep),
                 lower=c(rep(NA,Ji),yt_lower_dep), upper=c(rep(NA,Ji),yt_upper_dep))
p <- ggplot(df, aes(x, y)) + 
  geom_ribbon(aes(x, ymin=lower, ymax=upper), fill="gray", colour="gray") + 
  geom_line(size=2) + geom_hline(yintercept=16, colour='red', linetype='dashed') +
  theme_classic() + xlab("Visit") + ylab("Depression") + 
  theme(axis.title.x = element_text(size = 30)) + theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30))

name <- paste("real_state_pred_2.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

##########################################################################################
#                         Figure 8: Optimal cART Assignments                             #
##########################################################################################

# Individual S1

i <- 3 # subject index 
reward_coef <- c(0.4,0.5,0.1) # reward weights
uq_factor <- 0 # uncertainty quantification factors
source("SGD_Data.R") # Data processing for SGD 
sgd <- sgd_all$S1_uq0 

# (a) Optimal cART assignments for S1
dfa <- data.frame(visit = factor(rep(1:num_pred, 6), 
                                 levels=c("1", "2", "3", "4")),
                  Class = factor(rep(c("NRTI", "NNRTI", "PI", "INSTI", "EI", "Booster"), each=num_pred)),
                  num = c(c(2,2,2,2), c(0,0,0,0), c(0,0,0,0), c(1,1,1,1), c(0,1,1,1), c(0,0,0,0)))
dfb <- data.frame(drug_name = c(c("FTC","ABC","RAL"," "),c("FTC","ABC","RAL","MVC"),
                                c("FTC","ABC","RAL","MVC"),c("FTC","ABC","RAL","MVC")),
                  xpos = c(c(1, 1, 1, 1),c(2, 2, 2, 2),c(3, 3, 3, 3),c(4, 4, 4, 4)),
                  ypos = c(c(1, 2, 3, 4),c(1, 2, 3, 4),c(1, 2, 3, 4),c(1, 2, 3, 4)))
p <- ggplot(data=dfa, aes(x=visit, y=num, fill=Class)) + geom_bar(stat="identity") + scale_fill_viridis(discrete = TRUE) +
  geom_text(data=dfb, aes(y=ypos, x=xpos, label=drug_name), vjust=1.5, color="white", size=12, inherit.aes = FALSE) +
  theme_classic() + scale_y_continuous(limits=c(0,4)) + xlab("Visit") + ylab("Number of Drugs Used") + 
  scale_x_discrete(labels=c("1"="9","2"="10","3"="11","4"="12")) +
  theme(legend.title=element_text(size=30)) + theme(legend.text=element_text(size=30)) + 
  theme(axis.title.x = element_text(size=30)) + theme(axis.title.y = element_text(size=30)) +
  theme(axis.text.x = element_text(size=30)) + theme(axis.text.y = element_text(size=30)) + theme(legend.position='top')

name <- paste("real_optimal_1.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

# (c) Probability of selecting RAL for S1
k <- 5; iter <- which.max(unlist(sgd$R))
xx <- seq(0,60,by=0.1); yy <- seq(0,10,by=0.1)
zz <- matrix(NA,nrow=length(xx),ncol=length(yy))
for (i in 1:length(xx)){
  for (j in 1:length(yy)){
    y_pred_ij <- c(xx[i], yy[j], data$yp[3,Ji], data$yp[4,Ji])
    zz[i,j] <- softmax2(k, c(1,as.vector((y_pred_ij-y_mean)/y_std)), sgd$theta[[iter]]$theta2[[k]])[3] # INSTI: RAL
  }
}

df <- expand.grid(x=xx, y=yy)
df$Pr <- as.vector(zz)
p <- ggplot(df, aes(x, y, z=Pr)) + geom_raster(data=df, aes(fill=Pr), show.legend = TRUE)+
  geom_contour(colour="white", size=1, binwidth=0.5) + scale_fill_continuous(type="viridis") + 
  theme_classic() + xlab("Depression") + ylab("Viral Load") + 
  theme(legend.title=element_text(size=30)) + theme(legend.text=element_text(size=30)) + 
  theme(axis.title.x = element_text(size=30)) + theme(axis.title.y = element_text(size=30)) +
  theme(axis.text.x = element_text(size=30)) + theme(axis.text.y = element_text(size=30)) + 
  theme(legend.key.height= unit(2, 'cm'))

name <- paste("real_contour_1.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

# Individual S2

i <- 84 # subject index 
reward_coef <- c(1/3,1/3,1/3) # reward weights
uq_factor <- 0 # uncertainty quantification factor
source("SGD_Data.R") # Data processing for SGD 
sgd <- sgd_all$S2_uq0

# (b) Optimal cART assignments for S2
dfa <- data.frame(visit = factor(rep(1:num_pred, 6), 
                                 levels=c("1", "2", "3", "4")),
                  Class = rep(c("NRTI", "NNRTI", "PI", "INSTI", "EI", "Booster"), each=num_pred),
                  num = c(c(2,2,2,2), c(0,0,0,0), c(1,1,1,1), c(1,1,1,1), c(0,0,0,0), c(1,2,2,2)))
dfb <- data.frame(drug_name = c(c("ATV","FTC","ABC","DTG","RTV"," "),c("ATV","FTC","ABC","EVG","RTV","COBI"),
                                c("ATV","FTC","ABC","EVG","RTV","COBI"),c("LPV","FTC","ABC","EVG","RTV","COBI")),
                  xpos = c(c(1, 1, 1, 1, 1, 1),c(2, 2, 2, 2, 2, 2),c(3, 3, 3, 3, 3, 3),c(4, 4, 4, 4, 4, 4)),
                  ypos = c(c(1, 2, 3, 4, 5, 6),c(1, 2, 3, 4, 5, 6),c(1, 2, 3, 4, 5, 6),c(1, 2, 3, 4, 5, 6)))
p <- ggplot(data=dfa, aes(x=visit, y=num, fill=Class)) + geom_bar(stat="identity") + scale_fill_viridis(discrete = TRUE) + 
  geom_text(data=dfb, aes(y=ypos, x=xpos, label=drug_name), vjust=1.5, color="white", size=12, inherit.aes = FALSE) +
  theme_classic() + scale_y_continuous(limits=c(0,6)) + xlab("Visit") + ylab("Number of Drugs Used") + 
  scale_x_discrete(labels=c("1"="32","2"="33","3"="34","4"="35")) + 
  theme(legend.title=element_text(size=30)) + theme(legend.text=element_text(size=30)) + 
  theme(axis.title.x = element_text(size=30)) + theme(axis.title.y = element_text(size=30)) +
  theme(axis.text.x = element_text(size=30)) + theme(axis.text.y = element_text(size=30)) + theme(legend.position='top')

name <- paste("real_optimal_2.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

# (d) Probability of selecting LPV for S2
k <- 4; iter <- which.max(unlist(sgd$R))
xx <- seq(0,60,by=0.1); yy <- seq(0,10,by=0.1)
zz <- matrix(NA,nrow=length(xx),ncol=length(yy))
for (i in 1:length(xx)){
  for (j in 1:length(yy)){
    y_pred_ij <- c(xx[i], yy[j], data$yp[3,Ji], data$yp[4,Ji])
    zz[i,j] <- softmax2(k, c(1,as.vector((y_pred_ij-y_mean)/y_std)), sgd$theta[[iter]]$theta2[[k]])[3] # PI: LPV
  }
}

df <- expand.grid(x=xx, y=yy)
df$Pr <- as.vector(zz)
p <- ggplot(df, aes(x, y, z=Pr)) + geom_raster(data=df, aes(fill=Pr), show.legend = TRUE)+
  geom_contour(colour="white", size=1, binwidth=0.5) + scale_fill_continuous(type="viridis") +
  theme_classic() + xlab("Depression") + ylab("Viral Load") + 
  theme(legend.title=element_text(size=30)) + theme(legend.text=element_text(size=30)) + 
  theme(axis.title.x = element_text(size=30)) + theme(axis.title.y = element_text(size=30)) +
  theme(axis.text.x = element_text(size=30)) + theme(axis.text.y = element_text(size=30)) + 
  theme(legend.key.height= unit(2, 'cm'))

name <- paste("real_contour_2.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

##########################################################################################
#             Figure 9: SGD Results (uncertainty-penalized policy optimization)          #
##########################################################################################

# Individual P1

i <- 4 # subject index 
reward_coef <- c(0.1,0.1,0.8) # reward weights
uq_factor <- c(0,0.05,0.1) # uncertainty quantification factor
source("SGD_Data.R") # Data processing for SGD 

# (a) Observed eGFR for P1
df <- data.frame(x=1:Ji, y=data$yp[renal_index,])
p <- ggplot(df, aes(x, y)) + geom_line(size=2) + geom_hline(yintercept=60, colour='red', linetype='dashed') +
  theme_classic() + xlab("Visit") + ylab("eGFR") + 
  theme(axis.title.x = element_text(size = 30)) + theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30))

name <- paste("real_uq_1.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

# (b) Predicted eGFR for P1
y_pred <- list(); Z_pred <- list()
yt <- array(NA, dim=c(num_samples, num_pred, M))
y_pred2 <- list(); Z_pred2 <- list()
yt2 <- array(NA, dim=c(num_samples, num_pred, M))
y_pred3 <- list(); Z_pred3 <- list()
yt3 <- array(NA, dim=c(num_samples, num_pred, M))
# generate data for each posterior samples
for (nit in 1:num_samples){
  
  # lambda = 0
  data_pred <- data
  data_pred$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred[[nit]] <- list(); Z_pred[[nit]] <- list()
  Z_pred[[nit]][[1]] <- c("ABC","EFV","FTC")
  Z_pred[[nit]][[2]] <- c("ABC","EFV","FTC")
  Z_pred[[nit]][[3]] <- c("ABC","EFV","FTC")
  Z_pred[[nit]][[4]] <- c("ABC","EFV","FTC")
  
  # lambda = 0.05
  data_pred2 <- data
  data_pred2$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred2[[nit]] <- list(); Z_pred2[[nit]] <- list()
  Z_pred2[[nit]][[1]] <- c("3TC","ABC","EFV")
  Z_pred2[[nit]][[2]] <- c("3TC","ABC","EFV")
  Z_pred2[[nit]][[3]] <- c("3TC","ABC","EFV")
  Z_pred2[[nit]][[4]] <- c("3TC","ABC","EFV")
  
  # lambda = 0.1
  data_pred3 <- data
  data_pred3$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred3[[nit]] <- list(); Z_pred3[[nit]] <- list()
  Z_pred3[[nit]][[1]] <- c("3TC","ABC","ATV","RTV")
  Z_pred3[[nit]][[2]] <- c("3TC","ABC","ATV","RTV")
  Z_pred3[[nit]][[3]] <- c("3TC","ABC","ATV","RTV")
  Z_pred3[[nit]][[4]] <- c("3TC","ABC","ATV","RTV")
  
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
    # generate y_{Ji+1} conditional on the current regimen Z_{Ji+1}
    state_list <- generate_state_pred(data_pred, Ji+j-1, Z_pred[[nit]][[j]], tn[j], phi_nit)
    y_pred[[nit]][[j]] <- state_list$yn; data_pred <- state_list$data; yt[nit,j,] <- y_pred[[nit]][[j]]
    state_list <- generate_state_pred(data_pred2, Ji+j-1, Z_pred2[[nit]][[j]], tn[j], phi_nit)
    y_pred2[[nit]][[j]] <- state_list$yn; data_pred2 <- state_list$data; yt2[nit,j,] <- y_pred2[[nit]][[j]]
    state_list <- generate_state_pred(data_pred3, Ji+j-1, Z_pred3[[nit]][[j]], tn[j], phi_nit)
    y_pred3[[nit]][[j]] <- state_list$yn; data_pred3 <- state_list$data; yt3[nit,j,] <- y_pred3[[nit]][[j]]
  }
}

df <- data.frame(x=factor(c(rep(22,num_samples), rep(23,num_samples), rep(24,num_samples), rep(25, num_samples),
                            rep(22,num_samples), rep(23,num_samples), rep(24,num_samples), rep(25, num_samples),
                            rep(22,num_samples), rep(23,num_samples), rep(24,num_samples), rep(25, num_samples)),levels=c(22,23,24,25)), 
                 y=c(yt[,1,3], yt[,2,3], yt[,3,3], yt[,4,3], yt2[,1,3], yt2[,2,3], yt2[,3,3], yt2[,4,3],yt3[,1,3], yt3[,2,3], yt3[,3,3], yt3[,4,3]),
                 z=factor(c(rep(c("FTC+ABC+EFV"),num_samples*4), rep(c("3TC+ABC+EFV"),num_samples*4),rep(c("3TC+ABC+ATV+RTV"),num_samples*4)), 
                          levels=c("FTC+ABC+EFV","3TC+ABC+EFV","3TC+ABC+ATV+RTV")))
p <- ggplot(df, aes(x, y, fill=z)) + 
  stat_boxplot(geom ='errorbar', width = 0.25, position=position_dodge(0.75)) + 
  geom_boxplot(outlier.shape=NA,na.rm=TRUE) + theme_classic() + xlab("Visit") + ylab("eGFR") +
  geom_hline(yintercept=60, colour='red', linetype='dashed') + 
  scale_fill_grey(start = 1, end = 0.35, name=expression(lambda), labels=c(0,0.05,0.1)) +
  theme(axis.title.x = element_text(size = 30)) + theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30)) +
  theme(legend.title=element_text(size=30), legend.text=element_text(size=25)) 

name <- paste("real_uq_2.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()
