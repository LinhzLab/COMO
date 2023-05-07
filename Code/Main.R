rm(list=ls())
library(sp)
library(spdep)
library(splm)
library(Matrix)
library(nloptr)
####Generate simulation data###
K = 10;n <- 49*K;n1 <- 1000;N <- 1000
rho <- -0.9;beta <- c(0.3,-0.3,1);rho_hat <- rep(0,N)
beta_hat <- matrix(rep(0, 3*N), ncol=3,byrow=T);response <- matrix(0,n,1)
I_2 <- `diag<-`(matrix(0, K, K), 1)
I <- `diag<-`(matrix(0, n, n), 1)
Q = rep(1,n)%*%t(rep(1,n))-I
require(maptools)
col.gal.nb <- read.gal(system.file("etc/weights/columbus.gal",
                                   package="spdep")[1])
w<- nb2mat(col.gal.nb, zero.policy=TRUE)
w <- kronecker(I_2,w)
w_nb <- mat2listw (w)$neighbours
#################################
##Conduct simulation##
source("Como.R")
for(m in 1:n1){
  set.seed(m)
  epsilon <- rnorm(n,0,0.5)
  epsilon <- as.matrix(epsilon)
  x <- cbind(rnorm(n,2,0.25),rnorm(n,2,0.25),rep(1,n))
  y = solve(I-rho*w)%*%x%*%beta +solve(I-w*rho)%*%epsilon 
  data<-data.frame(y=y,x=x)
  wy =  w%*%y 
  estimation <- stsls(y~x[,1:2],data,listw = nb2listw(w_nb))
  rho_inital<- estimation$coefficients[1]
  beta_inital<- c(estimation$coefficients[3:4],estimation$coefficients[2])
  res <- COMO_algo(x,y,wy,Q,beta_inital,rho_inital)
  rho_hat[m] <-  res$rho_1
  beta_hat[m,1:3] <- res$solution_2
  print(m)
}

####Summarize results####
#hist(rho_hat)
mes_beta <- function(x){
  return(mean(x*x))
}
index = which(rho_hat< -0&rho_hat> - 2)
beta_true <-matrix(rep(beta,1),ncol = 3,byrow=T)
Mean_rho <- mean(rho_hat[index]) 
SD_rho <- sd(rho_hat[index])
MSE_rho <-SD_rho^2+ (Mean_rho-rho)^2

Mean_beta <- apply(beta_hat[index,],2,mean)
SD_beta <- matrix(rep(0,3),ncol = 3,byrow=T)
SD_beta[,1] <- sd(beta_hat[index,1])
SD_beta[,2] <- sd(beta_hat[index,2])
SD_beta[,3] <- sd(beta_hat[index,3]) 

bias_beta <- apply(Mean_beta-beta_true,2,mes_beta)
MSE_beta <-  apply(SD_beta,2,mes_beta)+bias_beta
len <- length(rho_hat[rho_hat!=0])
data <- data.frame(bias_rho = abs(Mean_rho-rho),SD_rho,sqrt(MSE_rho),bias_beta1 = sqrt(bias_beta[1]),sd_beta1 = SD_beta[1],rmse_beta1 = sqrt(MSE_beta[1]),bais_beta2 = sqrt(bias_beta[2]),sd_beta2 = SD_beta[2],rmse_beta2 = sqrt(MSE_beta[2]),bias_beta3 = sqrt(bias_beta[3]),sd_beta3 = SD_beta[3],rmse_beta3 = sqrt(MSE_beta[3]))
data <- round(data,digits=3)  
print(data)
path <- paste0("COMO_n_",n,".Rdata",sep = '')
save.image(file=path)


