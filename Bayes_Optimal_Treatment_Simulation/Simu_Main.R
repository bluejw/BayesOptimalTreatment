##########################################################################################
#                                   Simulation Study                                     #
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
load("Simu_Data_Preprocess.Rdata") # load the preprocessed data 
source("Simu_Data_Generate.R") # generate simulation truths 

##########################################################################################
#                           Run MCMC on Simulated Data                                   #
##########################################################################################

source("MCMC_Main.R") # MCMC main function

Nit <- 10000 # number of total MCMC intertaions 
burn.in <- 5000 # the burn-in iterations
thin.fac <- 10 # thinning factor for post burn-in samples 
MCMC_Results <- MCMC_Main(Nit, burn.in, thin.fac) # about 40 mins 

##########################################################################################
#                           Run SGD on Simulated Data                                    #
##########################################################################################

source("SGD_R_Functions.R") # R functions for SGD
sourceCpp("SGD_Rcpp_Functions.cpp") # Rcpp functions for SGD
source("Drug_Similarity.R") # subset-tree drug similarity 
source("SGD_Main.R") # SGD main function
load("Simu_MCMC_Results.Rdata") # posterior samples 

i <- 1 # subject index 
phi <- post # posterior samples 
Niter <- 1000 # number of SGD iterations
num_samples <- 100 # number of posterior samples
num_pred <- 4 # number of forward predictions
step_size <- 0.1 # step size of gradient descent
reward_coef <- c(1/3,1/3,1/3) # reward weights
uq_factor <- 0 # uncertainty quantification factor
source("SGD_Data.R") # Data processing for SGD 
SGD_Results <- SGD_Main(i, phi, Niter, num_samples, num_pred, step_size, reward_coef, uq_factor) # about 1 hour

##########################################################################################
#                           MCMC Posterior Inference                                     #
##########################################################################################

mean <- NULL # posterior means of parameters 
mean$theta <- colMeans(post$theta)
mean$delta <- colMeans(post$delta)
mean$rho <- mean(post$rho)
mean$sigma2 <- colMeans(post$sigma2)

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
#                       Figure S1: Trace Plots for Parameters                            #
##########################################################################################

theme_update(plot.title = element_text(hjust = 0.5, size = 25))

df11 <- data.frame(iteration = seq(1, 500, by=1), value = post$theta[,1,1])
p11 <- ggplot(df11, aes(iteration, value)) + geom_line() +  
  geom_hline(aes(yintercept = beta[1,1]),  colour = "red") +
  labs(title=expression(beta['1,1'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 

df12 <- data.frame(iteration = seq(1, 500, by=1), value = post$theta[,2,3])
p12 <- ggplot(df12, aes(iteration, value)) + geom_line() +  
  geom_hline(aes(yintercept = beta[2,3]),  colour = "red") +
  labs(title=expression(beta['2,3'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 

df13 <- data.frame(iteration = seq(1, 500, by=1), value = post$theta[,3,2])
p13 <- ggplot(df13, aes(iteration, value)) + geom_line() +  
  geom_hline(aes(yintercept = beta[3,2]),  colour = "red") +
  labs(title=expression(beta['3,2'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 

df21 <- data.frame(iteration = seq(1, 500, by=1), value = post$theta[,1,S+2])
p21 <- ggplot(df21, aes(iteration, value)) + geom_line() +  
  geom_hline(aes(yintercept = gamma[1,2]),  colour = "red") +
  labs(title=expression(gamma['1,2'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 

df22 <- data.frame(iteration = seq(1, 500, by=1), value = post$theta[,2,S+3])
p22 <- ggplot(df22, aes(iteration, value)) + geom_line() +  
  geom_hline(aes(yintercept = gamma[2,3]),  colour = "red") +
  labs(title=expression(gamma['2,3'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 

df23 <- data.frame(iteration = seq(1, 500, by=1), value = post$theta[,3,S+4])
p23 <- ggplot(df23, aes(iteration, value)) + geom_line() +  
  geom_hline(aes(yintercept = gamma[3,4]),  colour = "red") +
  labs(title=expression(gamma['3,4'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 

df31 <- data.frame(iteration = seq(1, 500, by=1), value = post$delta[,2,1])
p31 <- ggplot(df31, aes(iteration, value)) + geom_line() +  
  geom_hline(aes(yintercept = delta[2,1]),  colour = "red") +
  labs(title=expression(delta['2,1'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 

df32 <- data.frame(iteration = seq(1, 500, by=1), value = post$delta[,3,9])
p32 <- ggplot(df32, aes(iteration, value)) + geom_line() +  
  geom_hline(aes(yintercept = delta[3,9]),  colour = "red") +
  labs(title=expression(delta['3,9'])) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 

df33 <- data.frame(iteration = seq(1, 500, by=1), value = post$sigma2[,2])
p33 <- ggplot(df33, aes(iteration, value)) + geom_line() +  
  geom_hline(aes(yintercept = sigma2[2]),  colour = "red") +
  labs(title=expression(sigma[2]^{2})) + xlab("MCMC Iteration") + ylab("Value") + 
  theme(axis.title.x = element_text(size = 20)) + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) + theme(axis.text.y = element_text(size = 20)) 

name <- paste("trace_plot.pdf")
pdf(name, width = 18, height = 18, onefile = TRUE)
grid.arrange(p11, p12, p13, p21, p22, p23, p31, p32, p33, widths = c(1, 1, 1),
             layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9)))
dev.off()

##########################################################################################
#                  Figure S2: 95% Credible Intervals for Parameters                      #
##########################################################################################

df <- data.frame(par = rep(c(1,2,3,4,5),3), value = as.vector(c(t(theta[,c(2,3,4,5)]), delta[,4])),
                 lower = as.vector(c(t(ci$theta[,c(2,3,4,5),1]), ci$delta[,4,1])), 
                 upper = as.vector(c(t(ci$theta[,c(2,3,4,5),2]), ci$delta[,4,2])),
                 state = factor(c(rep("1",5), rep("2",5), rep("3",5)), levels=c("1","2","3")))
p <- ggplot(df, aes(par, value, colour = state)) + coord_flip() + 
  theme_classic() + geom_hline(yintercept=0, colour="darkgray") + 
  geom_pointrange(aes(ymin = lower, ymax = upper), position=position_dodge(0.5), size=1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0, position=position_dodge(0.5), size=1.5) + 
  labs(caption = " ", x = NULL, y = NULL, col = "m") + 
  scale_x_discrete(limit=c(1:5), labels = c(expression(beta['m,2']), expression(beta['m,3']), 
    expression(gamma['m,1']), expression(gamma['m,2']), expression(delta['m,4']))) +
  theme(legend.title=element_text(size=30)) + theme(legend.text=element_text(size=30)) + 
  theme(axis.title.x = element_text(size=30)) + theme(axis.title.y = element_text(size=30)) +
  theme(axis.text.x = element_text(size=30)) + theme(axis.text.y = element_text(size=30))

name <- paste("simu_esti.pdf")
pdf(name, width = 18, height = 8, onefile = TRUE)
p
dev.off()

##########################################################################################
#                        Table S3: Mean Squared Error                                    #
##########################################################################################

mse <- NULL; sd <- NULL
# theta
theta_array <- array(rep(theta,500),dim=c(M,S_tilde,500))
theta_array <- aperm(theta_array,c(3,1,2))
mse$theta <- colMeans((post$theta - theta_array)^2,dims=1)
sd$theta <- apply((post$theta - theta_array)^2, c(2,3), sd)
# delta
delta_array <- array(rep(delta,500),dim=c(M,Nk_sum,500))
delta_array <- aperm(delta_array,c(3,1,2))
mse$delta <- colMeans((post$delta - delta_array)^2,dims=1)
sd$delta <- apply((post$delta - delta_array)^2, c(2,3), sd)
# rho 
mse$rho <- mean((post$rho - rho)^2)
sd$rho <- sd((post$rho - rho)^2)
# sigma2
sigma2_mat <- array(rep(sigma2,500),dim=c(M,500))
sigma2_mat <- t(sigma2_mat)
mse$sigma2 <- colMeans((post$sigma2 - sigma2_mat)^2,dims=1)
sd$sigma2 <- apply((post$sigma2 - sigma2_mat)^2, 2, sd)

##########################################################################################
#                   Figure 3: SGD Results (uncertainty is not penalized)                 #
##########################################################################################

load("Simu_SGD_Results.Rdata") # SGD results for simulation study

# Individual I1

i <- 1 # subject index 
reward_coef <- c(1/3,1/3,1/3) # reward weights
uq_factor <- 0 # uncertainty quantification factor
source("SGD_Data.R") # Data processing for SGD 
sgd <- sgd_all$I1_uq0

# (a) Expected reward for I1
df <- data.frame(iteration=seq(1,1000,by=1), reward=unlist(sgd$R))
p <- ggplot(df, aes(iteration, reward)) + geom_line(colour = "black") +  
  xlab("SGD Iteration") + ylab("Expected Reward") + 
  theme_classic() + xlab("SGD Iteration") + ylab("Expected Reward") + 
  theme(axis.title.x = element_text(size = 30)) + theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30)) 

name <- paste("simu_sgd_reward_1.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

# (c) Viral load for I1
y_pred <- list(); Z_pred <- list()
yt <- array(NA, dim=c(num_samples, num_pred, M))
# generate data for each posterior samples
for (nit in 1:num_samples){
  
  data_pred <- data
  data_pred$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred[[nit]] <- list(); Z_pred[[nit]] <- list()
  Z_pred[[nit]][[1]] <- c("ABC","FTC","MVC","NVP")
  Z_pred[[nit]][[2]] <- c("ABC","FTC","MVC","NVP")
  Z_pred[[nit]][[3]] <- c("ABC","FTC","MVC","NVP")
  Z_pred[[nit]][[4]] <- c("ABC","ETR","FTC","MVC")
  
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
  scale_x_continuous(labels = as.character(1:length(df$x)), breaks = 1:length(df$x)) + 
  scale_y_continuous(limits=c(-2.75,7.5))

name <- paste("simu_state_pred_1.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

# Individual I2

i <- 2 # subject index 
reward_coef <- c(1/3,1/3,1/3) # reward weights
uq_factor <- 0 # uncertainty quantification factor
source("SGD_Data.R") # Data processing for SGD 
sgd <- sgd_all$I2_uq0

# (b) Expected reward for I2
df <- data.frame(iteration=seq(1,1000,by=1), reward=unlist(sgd$R))
p <- ggplot(df, aes(iteration, reward)) + geom_line(colour = "black") +  
  xlab("SGD Iteration") + ylab("Expected Reward") +
  theme_classic() + xlab("SGD Iteration") + ylab("Expected Reward") + 
  theme(axis.title.x = element_text(size = 30)) + theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30)) 

name <- paste("simu_sgd_reward_2.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

# (d) eGFR for I2
y_pred <- list(); Z_pred <- list()
yt <- array(NA, dim=c(num_samples, num_pred, M))
# generate data for each posterior samples
for (nit in 1:num_samples){
  
  data_pred <- data
  data_pred$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred[[nit]] <- list(); Z_pred[[nit]] <- list()
  Z_pred[[nit]][[1]] <- c("ABC","FTC","RPV")
  Z_pred[[nit]][[2]] <- c("ABC","FTC","RPV")
  Z_pred[[nit]][[3]] <- c("ABC","FTC","RPV")
  Z_pred[[nit]][[4]] <- c("ABC","FTC","RPV")
  
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

y_renal <- data$yp[renal_index,]
yt_mean_renal<- colMeans(yt[,,renal_index])
yt_lower_renal <- apply(yt[,,renal_index], 2, function(x){quantile(x, 0.025)})
yt_upper_renal <- apply(yt[,,renal_index], 2, function(x){quantile(x, 0.975)})

df <- data.frame(x=seq(1,Ji+num_pred,by=1), y=c(y_renal, yt_mean_renal),
                 lower=c(rep(NA,Ji),yt_lower_renal), upper=c(rep(NA,Ji),yt_upper_renal))
p <- ggplot(df, aes(x, y)) + 
  geom_ribbon(aes(x, ymin=lower, ymax=upper), fill="gray", colour="gray") + 
  geom_line(size=2) + geom_hline(yintercept=60, colour='red', linetype='dashed') +
  theme_classic() + xlab("Visit") + ylab("eGFR") + 
  theme(axis.title.x = element_text(size = 30)) + theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30))

name <- paste("simu_state_pred_2.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

##########################################################################################
#             Figure 4: SGD Results (uncertainty-penalized policy optimization)          #
##########################################################################################

# Individual I1

i <- 1 # subject index 
reward_coef <- c(1/3,1/3,1/3) # reward weights
uq_factor <- c(0,0.1,0.25,0.5) # uncertainty quantification factors
source("SGD_Data.R") # Data processing for SGD 

# (a) Predicted depression score for I1
y_pred <- list(); Z_pred <- list()
yt <- array(NA, dim=c(num_samples, num_pred, M))
y_pred2 <- list(); Z_pred2 <- list()
yt2 <- array(NA, dim=c(num_samples, num_pred, M))
y_pred3 <- list(); Z_pred3 <- list()
yt3 <- array(NA, dim=c(num_samples, num_pred, M))
y_pred4 <- list(); Z_pred4 <- list()
yt4 <- array(NA, dim=c(num_samples, num_pred, M))
# generate data for each posterior samples
for (nit in 1:num_samples){
  
  data_pred <- data
  data_pred$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred[[nit]] <- list(); Z_pred[[nit]] <- list()
  Z_pred[[nit]][[1]] <- c("ABC","FTC","MVC","NVP")
  Z_pred[[nit]][[2]] <- c("ABC","FTC","MVC","NVP")
  Z_pred[[nit]][[3]] <- c("ABC","FTC","MVC","NVP")
  Z_pred[[nit]][[4]] <- c("ABC","ETR","FTC","MVC")
  
  data_pred2 <- data
  data_pred2$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred2[[nit]] <- list(); Z_pred2[[nit]] <- list()
  Z_pred2[[nit]][[1]] <- c("3TC","ABC","NVP")
  Z_pred2[[nit]][[2]] <- c("3TC","ABC","NVP")
  Z_pred2[[nit]][[3]] <- c("3TC","ABC","NVP")
  Z_pred2[[nit]][[4]] <- c("3TC","ABC","NVP")
  
  data_pred3 <- data
  data_pred3$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred3[[nit]] <- list(); Z_pred3[[nit]] <- list()
  Z_pred3[[nit]][[1]] <- c("3TC","ABC","ATV","NVP","RTV")
  Z_pred3[[nit]][[2]] <- c("3TC","ABC","ATV","NVP","RTV")
  Z_pred3[[nit]][[3]] <- c("3TC","ABC","ATV","NVP","RTV")
  Z_pred3[[nit]][[4]] <- c("3TC","ABC","ATV","NVP","RTV")
  
  data_pred4 <- data
  data_pred4$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred4[[nit]] <- list(); Z_pred4[[nit]] <- list()
  Z_pred4[[nit]][[1]] <- c("3TC","ABC","ATV","NVP","RTV")
  Z_pred4[[nit]][[2]] <- c("3TC","ABC","ATV","NVP","RTV")
  Z_pred4[[nit]][[3]] <- c("3TC","ABC","ATV","NVP","RTV")
  Z_pred4[[nit]][[4]] <- c("3TC","ABC","ATV","NVP","RTV")
  
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
    state_list <- generate_state_pred(data_pred4, Ji+j-1, Z_pred4[[nit]][[j]], tn[j], phi_nit)
    y_pred4[[nit]][[j]] <- state_list$yn; data_pred4 <- state_list$data; yt4[nit,j,] <- y_pred4[[nit]][[j]]
  }
}

index <- dep_index
df <- data.frame(x=factor(c(rep(8,num_samples), rep(9,num_samples), rep(10,num_samples), rep(11,num_samples),
                            rep(8,num_samples), rep(9,num_samples), rep(10,num_samples), rep(11,num_samples),
                            rep(8,num_samples), rep(9,num_samples), rep(10,num_samples), rep(11,num_samples),
                            rep(8,num_samples), rep(9,num_samples), rep(10,num_samples), rep(11,num_samples)),levels=c(8,9,10,11)), 
                 y=c(yt[,1,index], yt[,2,index], yt[,3,index], yt[,4,index], yt2[,1,index], yt2[,2,index], yt2[,3,index], yt2[,4,index],
                     yt3[,1,index], yt3[,2,index], yt3[,3,index], yt3[,4,index], yt4[,1,index], yt4[,2,index], yt4[,3,index], yt4[,4,index]),
                 z=c(rep(c("1"),num_samples*4), rep(c("2"),num_samples*4), rep(c("3"),num_samples*4), rep(c("4"),num_samples*4)))
p <- ggplot(df, aes(x, y, fill=z)) + 
  stat_boxplot(geom ='errorbar', width = 0.25, position=position_dodge(0.75)) + 
  geom_boxplot(outlier.shape=NA,na.rm=TRUE) + theme_classic() + xlab("Visit") + ylab("Depression") + 
  scale_fill_grey(start = 1, end = 0.35, name=expression(lambda), labels=c(0, 0.1, 0.25, 0.5)) +
  theme(axis.title.x = element_text(size = 30)) + theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30)) +
  theme(legend.title=element_text(size=30), legend.text=element_text(size=25)) +
  geom_hline(yintercept=16, colour='red', linetype='dashed')

name <- paste("simu_uq_1.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

# Individual I2

i <- 2 # subject index 
reward_coef <- c(1/3,1/3,1/3) # reward weights
uq_factor <- c(0,0.1,0.25,0.5) # uncertainty quantification factors
source("SGD_Data.R") # Data processing for SGD 

# (b) Predicted depression score for I2
y_pred <- list(); Z_pred <- list()
yt <- array(NA, dim=c(num_samples, num_pred, M))
y_pred2 <- list(); Z_pred2 <- list()
yt2 <- array(NA, dim=c(num_samples, num_pred, M))
y_pred3 <- list(); Z_pred3 <- list()
yt3 <- array(NA, dim=c(num_samples, num_pred, M))
y_pred4 <- list(); Z_pred4 <- list()
yt4 <- array(NA, dim=c(num_samples, num_pred, M))
# generate data for each posterior samples
for (nit in 1:num_samples){
  
  data_pred <- data
  data_pred$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred[[nit]] <- list(); Z_pred[[nit]] <- list()
  Z_pred[[nit]][[1]] <- c("ABC","FTC","RPV")
  Z_pred[[nit]][[2]] <- c("ABC","FTC","RPV")
  Z_pred[[nit]][[3]] <- c("ABC","FTC","RPV")
  Z_pred[[nit]][[4]] <- c("ABC","FTC","RPV")
  
  data_pred2 <- data
  data_pred2$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred2[[nit]] <- list(); Z_pred2[[nit]] <- list()
  Z_pred2[[nit]][[1]] <- c("ABC","FTC","RPV")
  Z_pred2[[nit]][[2]] <- c("ABC","FTC","RPV")
  Z_pred2[[nit]][[3]] <- c("ABC","FTC","RPV")
  Z_pred2[[nit]][[4]] <- c("ABC","FTC","RPV")
  
  data_pred3 <- data
  data_pred3$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred3[[nit]] <- list(); Z_pred3[[nit]] <- list()
  Z_pred3[[nit]][[1]] <- c("EFV","FTC","TDF")
  Z_pred3[[nit]][[2]] <- c("EFV","FTC","TDF")
  Z_pred3[[nit]][[3]] <- c("EFV","FTC","TDF")
  Z_pred3[[nit]][[4]] <- c("EFV","FTC","TDF")
  
  data_pred4 <- data
  data_pred4$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred4[[nit]] <- list(); Z_pred4[[nit]] <- list()
  Z_pred4[[nit]][[1]] <- c("EFV","FTC","TDF")
  Z_pred4[[nit]][[2]] <- c("EFV","FTC","TDF")
  Z_pred4[[nit]][[3]] <- c("EFV","FTC","TDF")
  Z_pred4[[nit]][[4]] <- c("EFV","FTC","TDF")
  
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
    state_list <- generate_state_pred(data_pred4, Ji+j-1, Z_pred4[[nit]][[j]], tn[j], phi_nit)
    y_pred4[[nit]][[j]] <- state_list$yn; data_pred4 <- state_list$data; yt4[nit,j,] <- y_pred4[[nit]][[j]]
  }
}

index <- dep_index
df <- data.frame(x=factor(c(rep(22,num_samples), rep(23,num_samples), rep(24,num_samples), rep(25,num_samples),
                            rep(22,num_samples), rep(23,num_samples), rep(24,num_samples), rep(25,num_samples),
                            rep(22,num_samples), rep(23,num_samples), rep(24,num_samples), rep(25,num_samples),
                            rep(22,num_samples), rep(23,num_samples), rep(24,num_samples), rep(25,num_samples)),levels=c(22,23,24,25)), 
                 y=c(yt[,1,index], yt[,2,index], yt[,3,index], yt[,4,index], yt2[,1,index], yt2[,2,index], yt2[,3,index], yt2[,4,index],
                     yt3[,1,index], yt3[,2,index], yt3[,3,index], yt3[,4,index], yt4[,1,index], yt4[,2,index], yt4[,3,index], yt4[,4,index]),
                 z=c(rep(c("1"),num_samples*4), rep(c("2"),num_samples*4), rep(c("3"),num_samples*4), rep(c("4"),num_samples*4)))
p <- ggplot(df, aes(x, y, fill=z)) + 
  stat_boxplot(geom ='errorbar', width = 0.25, position=position_dodge(0.75)) + 
  geom_boxplot(outlier.shape=NA,na.rm=TRUE) + theme_classic() + xlab("Visit") + ylab("Depression") + 
  scale_fill_grey(start = 1, end = 0.35, name=expression(lambda), labels=c(0, 0.1, 0.25, 0.5)) +
  theme(axis.title.x = element_text(size = 30)) + theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 30)) + theme(axis.text.y = element_text(size = 30)) +
  theme(legend.title=element_text(size=30), legend.text=element_text(size=25)) +
  geom_hline(yintercept=16, colour='red', linetype='dashed')

name <- paste("simu_uq_2.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

##########################################################################################
#                             Figure 5: Model Comparison                                 #
##########################################################################################

# Individual I1

i <- 1 # subject index 
reward_coef <- c(1/3,1/3,1/3) # reward weights
source("SGD_Data.R") # Data processing for SGD 

# state prediction
y_pred <- list(); Z_pred <- list()
yt <- array(NA, dim=c(num_samples, num_pred, M))
y_pred2 <- list(); Z_pred2 <- list()
yt2 <- array(NA, dim=c(num_samples, num_pred, M))
y_pred3 <- list(); Z_pred3 <- list()
yt3 <- array(NA, dim=c(num_samples, num_pred, M))
# generate data for each posterior samples
for (nit in 1:num_samples){
  
  # proposed method
  data_pred <- data 
  data_pred$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred[[nit]] <- list(); Z_pred[[nit]] <- list()
  Z_pred[[nit]][[1]] <- c("ABC","FTC","MVC","NVP")
  Z_pred[[nit]][[2]] <- c("ABC","FTC","MVC","NVP")
  Z_pred[[nit]][[3]] <- c("ABC","FTC","MVC","NVP")
  Z_pred[[nit]][[4]] <- c("ABC","ETR","FTC","MVC")
  
  # no switch method 
  data_pred2 <- data
  data_pred2$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred2[[nit]] <- list(); Z_pred2[[nit]] <- list()
  Z_pred2[[nit]][[1]] <- c("3TC","AZT","NFV")
  Z_pred2[[nit]][[2]] <- c("3TC","AZT","NFV")
  Z_pred2[[nit]][[3]] <- c("3TC","AZT","NFV")
  Z_pred2[[nit]][[4]] <- c("3TC","AZT","NFV")
  
  # short term method
  data_pred3 <- data
  data_pred3$yp <- phi$y[nit,i,1:M,1:Ji]
  y_pred3[[nit]] <- list(); Z_pred3[[nit]] <- list()
  Z_pred3[[nit]][[1]] <- c("3TC","ABC","ETR")
  Z_pred3[[nit]][[2]] <- c("3TC","ABC","ETR")
  Z_pred3[[nit]][[3]] <- c("3TC","ABC","ETR")
  Z_pred3[[nit]][[4]] <- c("3TC","ABC","ETR")
  
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

# (a) Predicted depression score for I1
df <- data.frame(x=factor(c(rep(8,num_samples), rep(9,num_samples), rep(10,num_samples), rep(11,num_samples),
                            rep(8,num_samples), rep(9,num_samples), rep(10,num_samples), rep(11,num_samples),
                            rep(8,num_samples), rep(9,num_samples), rep(10,num_samples), rep(11,num_samples)),levels=c(8,9,10,11)), 
                 y=c(yt[,1,1], yt[,2,1], yt[,3,1], yt[,4,1], 
                     yt2[,1,1], yt2[,2,1], yt2[,3,1], yt2[,4,1], 
                     yt3[,1,1], yt3[,2,1], yt3[,3,1], yt3[,4,1]),
                 z=factor(c(rep(c("Proposed"),num_samples*4), rep(c("No Switch"),num_samples*4), rep(c("Short Term"),num_samples*4)),
                          level=c("No Switch","Short Term","Proposed")))
p <- ggplot(df, aes(x, y, fill=z)) + 
  stat_boxplot(geom ='errorbar', width = 0.25, position=position_dodge(0.75)) + 
  geom_boxplot(outlier.shape=NA,na.rm=TRUE) + theme_classic() + xlab("Visit") + ylab("Depression") + 
  scale_fill_manual(name="Method", values=c("white","grey70","grey40"), labels=c("No Switch","Short Term","Proposed")) +
  theme(axis.title.x = element_text(size = 30)) + theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 25)) + theme(axis.text.y = element_text(size = 25)) +
  theme(legend.title=element_text(size=25), legend.text=element_text(size=25)) + 
  geom_hline(yintercept=16, colour='red', linetype='dashed')

name <- paste("simu_comp_1.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()

# (b) Predicted viral load for I1
df <- data.frame(x=factor(c(rep(8,num_samples), rep(9,num_samples), rep(10,num_samples), rep(11,num_samples),
                            rep(8,num_samples), rep(9,num_samples), rep(10,num_samples), rep(11,num_samples),
                            rep(8,num_samples), rep(9,num_samples), rep(10,num_samples), rep(11,num_samples)),levels=c(8,9,10,11)), 
                 y=c(yt[,1,2], yt[,2,2], yt[,3,2], yt[,4,2], 
                     yt2[,1,2], yt2[,2,2], yt2[,3,2], yt2[,4,2], 
                     yt3[,1,2], yt3[,2,2], yt3[,3,2], yt3[,4,2]),
                 z=factor(c(rep(c("Proposed"),num_samples*4), rep(c("No Switch"),num_samples*4), rep(c("Short Term"),num_samples*4)),
                          level=c("No Switch","Short Term","Proposed")))
p <- ggplot(df, aes(x, y, fill=z)) + 
  stat_boxplot(geom ='errorbar', width = 0.25, position=position_dodge(0.75)) + 
  geom_boxplot(outlier.shape=NA,na.rm=TRUE) + theme_classic() + xlab("Visit") + ylab("Viral Load") + 
  scale_fill_manual(name="Method", values=c("white","grey70","grey40"), labels=c("No Switch","Short Term","Proposed")) +
  theme(axis.title.x = element_text(size = 30)) + theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 25)) + theme(axis.text.y = element_text(size = 25)) +
  theme(legend.title=element_text(size=25), legend.text=element_text(size=25)) +
  geom_hline(yintercept=log(20), colour='red', linetype='dashed')

name <- paste("simu_comp_2.pdf")
pdf(name, width = 10, height = 6, onefile = TRUE)
p
dev.off()
