################################################################################
#                                                                                     
#   Filename    :    simulation_study.R    												  
#                                                                                     
#   Project     :    BiomJ article "A marginalized zero-inflated negative binomial
#                    model for spatial data: modeling COVID-19 deaths in Georgia"   
#   Authors     :    Fedelis Mutiso ,  Hong Li ,  John L. Pearce ,  Sara E. Benjamin-Neelon ,  Noel T. Mueller ,  and Brian Neelon                                                              
#   Date        :    05.20.2024
#   Purpose     :    Reproduce Table 1 ,  Figure 3 ,  and Figure 4 of the manuscript
#																				  
#   R Version   :    4.2.2 (2022-10-31 ucrt)                                                                
#
#   Input data files  :    ordered_adjmat_ga.txt
#
#   Output data files :   Table_1.csv ,  Figure_3a.pdf ,  Figure_3b.pdf ,  Figures 4_ab.jpg ,  Figure 4_cd.jpg
#
################################################################################


##############################
# load the required packages #
##############################

library(BayesLogit)   
library(mvtnorm)	
library(MCMCpack)    
library(msm)          
library(spam)         
library(splines2)
library(fields)
library(tidyverse)
library(data.table)

library(tigris)
library(rgdal)
library(sf)
library(tmap)
library(RColorBrewer)

setwd('~/')                 # specify the full path to the location where the `bimj_files` folder containing the r scripts and data/results subfolders is stored



################################
# Adjacency Matrix Information #
################################

A <- matrix(scan("./bimj_files/data/ordered_adjmat_ga.txt") ,  159 ,  159)
m <- apply(A , 1 , sum)	                                    # No. neighbors for each county
ncounty <- n <- nrow(A)                                     # number of counties

adjid <- rep(1:ncounty ,  m)                                # Tells which elements of A belong to which block
adj <- rep(0 ,  sum(m))                                     # how many in total are adjacent
for (i in 1:ncounty) adj[adjid==i] <- which(A[i ,  ]==1)


set.seed(091322)

########################
# Temporal data        #
########################

ncounty <- n <- 159              # number of counties
lday <- 365
nis <- rep(lday ,  ncounty)      # number of observations for each county
id <- rep(1:ncounty ,  nis)      # county id indicator for each day
N <- length(id)
days <- rep(1:lday ,  ncounty) 
pop <- rep(sample(1500:1000000 ,  ncounty ,  replace = TRUE) ,  nis)
lpop <- log(pop)



#############################################################
# Fixed effect spline coefs for overall time trend
#############################################################

knots <- seq(14 ,  351 ,  14)
splines.tmp <- tmp <- bSpline(1:lday ,  knots = knots ,  intercept = T)
k <- ncol(tmp)
U <- apply(tmp ,  2 ,  rep ,  ncounty)

###########################
# spline coefficients     #
###########################

alpha1 <- 1*c(-1 ,  0.5 ,   1.62 ,   0.16 ,   1.51 ,   0.80 ,   0.73 ,  -0.68 ,  -0.42 ,  -1.06 ,  -0.26 ,  -0.80 ,  -0.68 ,  -0.91 ,  -0.87 ,  -1.92 ,  -1.79 ,  -1.32 ,  -0.51 , 
              0.48 ,   0.61 ,   1.13 ,   0.56 ,  0.42 ,  -0.23 ,  -0.15 ,  -1.05 ,   0.22 ,  0.5)


alpha2 <- 1*c(-1 ,  0.5 ,   1.62 ,   0.16 ,   1.51 ,   0.80 ,   0.73 ,  -0.68 ,  -0.42 ,  -1.06 ,  -0.26 ,  -0.80 ,  -0.68 ,  -0.91 ,  -0.87 ,  -1.92 ,  -1.79 ,  -1.32 ,  -0.51 , 
              0.48 ,   0.61 ,   1.13 ,   0.56 ,  0.42 ,  -0.23 ,  -0.15 ,  -1.05 ,   0.22 ,  0.5)


########################
#  Spatial Effects     #
########################
kappa <- .999999			  	                  # Spatial dependence parameter ~ 1 for intrinsic CAR
cov <- matrix(c(.5 ,  0.25 , 
              0.25 ,  0.75) , 2 , 2)                              # Conditional Cov of phi1 and phi2 given phi_{-i}

Q <- as.spam(diag(m))-kappa*as.spam(A)			     
covphi <- solve(Q)%x%cov			                 # covariance of phis (2n x 2n)
phi <- spam::rmvnorm(1 , sigma = covphi)		         # 1 x 2n matrix of spatial effects
phitmp <- matrix(phi ,  ncol = 2 ,  byrow = T)                   # n x 2 matrix of spatial effects

true.phi1 <- phi1 <- phitmp[ , 1]-mean(phitmp[ , 1])             # n x 1 phi1 vector -- Centered
true.phi2 <- phi2 <- phitmp[ , 2]-mean(phitmp[ , 2])             # n x 1 phi2 vector ,  etc.

Phi1 <- rep(phi1 , nis)
Phi2 <- rep(phi2 , nis)

true.phi <- cbind(true.phi1 , true.phi2)



###################
# binomial part   #
###################

# Fixed effects
x1 <- rep(rnorm(ncounty ,  mean = 0 ,  sd = 1) ,  nis)

beta1.true <- beta1 <- c(-0.25 ,  alpha1)                            # combine the spline and fixed effect coefficients together
X1 <- cbind(x1 ,  U)

X2 <- cbind(x1 ,  U)                                                 # same covariates for the mean and binary part
p <- ncol(X2)

eta1 <- X1%*%beta1+Phi1
psi <- c(exp(eta1)/(1+exp(eta1)))
l <- rbinom(N ,  1 ,  psi)                               # at risk indicator (= 1 if at risk )
N1 <- sum(l)                                             # total number of at risk subjects across counties
pstruct0 <- 1-mean(l)                                    # proportion of structural zeros in all counties



###################
# count part      #
###################

beta2.true <- beta2 <- c(-.5 ,  alpha2)                 # combine the spline and fixed effect coefficients together


nis1 <- tapply(l ,  id ,  sum)                          # the number of at risk observations in each county

eta2 <- X2[l==1 , ]%*%beta2+lpop[l==1]-log(10^6)+Phi2[l==1] 

r <- 1                                          # NB dispersion 
nu <- exp(eta2)                                 # marginal mean nu
mu <- nu/psi[l==1]                              # the mean of the positive count part



y <- rep(0 , N)                              # Response
y[l==1] <- rnbinom(N1 ,  r ,  mu = mu)       # If at risk ,  draw from NB
pzero <- length(y[y==0])/N                   # Proportion of zeros



##########   
# Priors #
##########

beta10 <- rep(0 , p)
beta20 <- rep(0 , p)
T0b1 <- diag(.01 , p)
T0b2 <- diag(.01 , p)         
s <- 0.002                                    # proposal variance



###################################
# Initialize different parameters #
###################################

beta1 <- rep(0 ,  p)
beta2 <- rep(0 ,  p)


phi_init <- spam::rmvnorm(1 ,  sigma = diag(.1 ,  2*ncounty))	    # Random effects
phi_init <- matrix(phi_init ,  ncol = 2 ,  byrow = T)               # n x 2 matrix of spatial effects
phi1 <- phi_init[ , 1]
phi2 <- phi_init[ , 2]

Phi1 <- rep(phi1 , nis)
Phi2 <- rep(phi2 , nis)

phibar1 <- phibar2 <- rep(0 ,  ncounty)      # mean of adjacent phi2's for updating phi2
for (j in 1:ncounty){
  phibar1[j] <- mean(phi1[adj[which(adjid==j)]])
  phibar2[j] <- mean(phi2[adj[which(adjid==j)]])
}

sphi1 <- sd(phi1)
tauphi1 <- 1/sphi1^2

sphi2 <- sd(phi2)
tauphi2 <- 1/sphi2^2

Sigmaphi <- cov(phi_init)
rho <- Sigmaphi[1 , 2]/sqrt(Sigmaphi[1 , 1]*Sigmaphi[2 , 2])

bcov <- diag(0.01 ,  p)
sb2 <- 0.04                  # tuning parameter for covariance of parameters for the mean
accb <- 0


s2 <- 0.17                   # tuning parameter for the variance of the random effects for the mean
u <- rep(0 ,  ncounty)       # keeps track of acceptance of phi2 for each county at each iteration
accu <- 0


r <- 1.2
accr <- 0  
y1 <- rep(0 , N)             # At risk indicator (this is W in paper)
y1[y>0] <- 1                 # If y>0 ,  then at risk w.p. 1
N0 <- length(y[y==0])        # Number of observed 0's
q <- rep(.5 , N)             # 1-p=1/(1+exp(X%*%alpha)) ,  used for updating y1


############
# Num Sims #
############
nsim <- 20000			          # Number of MCMC Iterations 
thin <- 10		                  # Thinning interval
burn <- 10000		                  # Burnin
lastit <- (nsim-burn)/thin	          # Last stored value

#########
# Store #
#########
Beta1 <- matrix(0 , lastit ,  p)
Beta2 <- matrix(0 , lastit ,  p)
BETA1 <- matrix(0 ,  nsim ,  p)
BETA2 <- matrix(0 ,  nsim ,  p)             # for updating beta2 proposal
R <- rep(0 , lastit)
Sigphi <- matrix(0 , lastit ,  4)
PHI1 <- PHI2 <- matrix(0 ,  lastit ,  ncounty)


########
# MCMC #
########

tmptime <- proc.time()

for (i in 1:nsim){
  
  #################
  # Update beta1  #
  #################
  
  eta1 <- X1%*%beta1+Phi1
  w <- rpg(N , 1 , eta1)
  z <- (y1-1/2)/w           
  vb <- solve(crossprod(sqrt(w)*X1)+T0b1)  
  mb <- vb%*%(T0b1%*%beta10+t(sqrt(w)*X1)%*%(sqrt(w)*(z-Phi1))) 
  beta1 <- c(spam::rmvnorm(1 , mb , vb))
  
  BETA1[i ,  ] <- beta1
  
  ################
  # update phi1  #
  ################
  
  priorprec <- 1/(sphi1^2*(1-rho^2))*Q                                             # Prior precision of phi1|phi2
  priormean <- rho*sphi1/sphi2*(phi2)                                              # Prior mean of phi1|phi2 --- just 0 since just one component and rho=0
  prec <- priorprec+as.spam(diag(tapply(w ,  id ,  sum) ,  ncounty ,  ncounty))
  mb <- c(priorprec%*%priormean)+tapply(w*(z-X1%*%beta1) , id , sum)               # note priormean is 0 and only data contributes to mb
  phi1 <- rmvnorm.canonical(1 ,  mb ,  prec)[1 , ]
  
  
  # center phi1 and update tauphi1

  phi1 <- phi1-mean(phi1)
  Phi1 <- rep(phi1 ,  nis)
  
  # update phibar1   (apply the sum to zero constraint)

  for (j in 1:ncounty) phibar1[j] <- mean(phi1[adj[which(adjid==j)]])
  
  
  
  ###################################
  # update the at-risk indicator    #
  ###################################
  
  
  eta1 <- X1%*%beta1+Phi1
  eta2 <- X2%*%beta2+lpop-log(10^6)+Phi2                         # Use all n observations
  psi <- c(pmax(0.01 , pmin(0.99 , exp(eta1)/(1+exp(eta1)))))    # at-risk probability
  nu = c(exp(eta2))
  mu <- nu/psi
  q <- r/(mu+r)                                    # Pr(y=0|y1=1)
  theta <- psi*(q^r)/(psi*(q^r)+1-psi)             # Conditional prob that y1=1 given y=0 -- i.e. Pr(chance zero|observed zero)
  y1[y==0] <- rbinom(N0 , 1 , theta[y==0])         # If y=0 ,  then draw a "chance zero" w.p. theta ,  otherwise y1=1
  N1 <- sum(y1)
  nis1 <- tapply(y1 , id , sum)
  
  ###############
  # Update r    #
  ###############
  
  rnew <- rtnorm(1 , r , sqrt(s) , lower = 0)       # Treat r as continuous
  ratio <- sum(dnbinom(y[y1==1] ,  rnew ,  mu = mu[y1==1] , log = T))-sum(dnbinom(y[y1==1] ,  r ,  mu = mu[y1==1] , log = T))+   
    dtnorm(r , rnew , sqrt(s) , 0 , log = T) - dtnorm(rnew , r , sqrt(s) , 0 , log = T)   # Uniform Prior for R 
  # Proposal not symmetric 
  if (log(runif(1))<ratio) {
    r <- rnew
    accr <- accr+1
  }
  
  
  ################
  # update beta2 #
  ################
  
  # Proposal
  beta2s <- c(beta2+rmvt(1 ,  df = 3 ,  sigma=sb2*bcov))
  eta2s <- X2[y1==1 , ]%*%beta2s+lpop[y1==1]-log(10^6)+Phi2[y1==1]
  nus <- c(exp(eta2s))
  mus <- nus/psi[y1==1]
  ls <- sum(dnbinom(y[y1==1] ,  r ,  mu = mus , log = T))                                   # Log like under proposal
  lps <- dmvnorm(beta2s ,  mean = rep(0 ,  p) ,  sigma = diag(100 ,  p) , log = T)          # Log prior proposal
  
  # Current
  eta2 <- X2[y1==1 , ]%*%beta2+lpop[y1==1]-log(10^6)+Phi2[y1==1]
  nu <- c(exp(eta2))
  mu <- nu/psi[y1==1]
  l <- sum(dnbinom(y[y1==1] ,  r ,  mu = mu , log = T))
  lp <- dmvnorm(beta2 ,  mean = rep(0 ,  p) ,  sigma=diag(100 ,  p) , log = T)
  
  ratio <- ls+lps-l-lp   # Note: Prior for r is lognormal ,  so log(r) is normal ,  
  
  if (log(runif(1))<ratio) {
    beta2 <- beta2s
    accb <- accb+1
  }
  
  BETA2[i , ] <- beta2
  if (i==nsim/2) bcov <- cov(BETA2[(nsim/4+1):(nsim/2) , ])                # Tune proposal covariance for Beta
  
  
  ###############
  # Update phi2 #
  ###############
  
  # Proposal
  
  phi2s <- phi2+s2*rt(ncounty ,  3)
  Phi2s <- rep(phi2s ,  nis)
  eta2s <- X2[y1==1 , ]%*%beta2+Phi2s[y1==1]+lpop[y1==1]-log(10^6)
  nus <- c(exp(eta2s))
  mus <- nus/psi[y1==1]
  Ls <- dnbinom(y[y1==1] ,  r ,  mu = mus ,  log = T)
  ls <- rep(0 ,  ncounty)  # account for empty counties - their likelihood contribution is 0
  ls[nis1>0] <- tapply(Ls ,  id[y1==1] ,  sum)
  
  # Current
  
  eta2 <- X2[y1==1 , ]%*%beta2+Phi2[y1==1]+lpop[y1==1]-log(10^6)
  nu <- c(exp(eta2))
  mu <- nu/psi[y1==1]
  L <- dnbinom(y[y1==1] ,  r ,  mu = mu ,  log = T)
  l <- rep(0 ,  ncounty)
  l[nis1>0] <- tapply(L ,  id[y1==1] ,  sum)
  
  # prior -- both current and proposal
  for (j in 1:ncounty){ 
    lp <- dnorm(phi2[j] ,  phibar2[j]+rho*sphi2/sphi1*(phi1[j]-phibar1[j]) ,  sqrt((sphi2^2*(1-rho^2)/m[j])) ,  log = T)
    lps <- dnorm(phi2s[j] ,  phibar2[j]+rho*sphi2/sphi1*(phi1[j]-phibar1[j]) ,  sqrt((sphi2^2*(1-rho^2)/m[j])) , log = T)
    ratio <- ls[j]+lps-l[j]-lp
    u[j] <- 1*(log(runif(1))<ratio)                                # keep track of acceptance of phi2 for each county
    if (log(runif(1))<ratio) {
      phi2[j] <- phi2s[j]
    }
    phibar2[j] <- mean(phi2[adj[which(adjid==j)]])                 # only updated if current phi2 is accepted
  }
  
  accu[i] <- mean(u)
  
  
  # Center phi2 
  
  phi2 <- phi2-mean(phi2) 
  Phi2 <- rep(phi2 ,  nis)
  
  #######################
  # Update Sigma.phi    #
  #######################
  
  phimat <- cbind(phi1 , phi2)
  Sigmaphi <- riwish(3+ncounty-1 , diag(2)+t(phimat)%*%Q%*%phimat)
  sphi1 <- sqrt(Sigmaphi[1 , 1])
  sphi2 <- sqrt(Sigmaphi[2 , 2])
  rho <- Sigmaphi[1 , 2]/sqrt(Sigmaphi[1 , 1]*Sigmaphi[2 , 2])
  
  # Store
  if (i> burn & i%%thin==0) {
    j <- (i-burn)/thin
    Beta1[j , ] <- beta1
    Beta2[j , ] <- beta2
    R[j] <- r		
    Sigphi[j , ] <- c(Sigmaphi)
    PHI1[j , ] <- phi1
    PHI2[j , ] <- phi2
  }
  
  if (i%%10==0){
    print(i)
    # print(c(beta1))
    # print(c(beta2))
    # print(c(Sigmaphi))
    # print(r)
  }
}

tot.time <- proc.time()-tmptime 


##################################
# posterior mean estimates       #
##################################

####################
# Binary component #
####################

beta1.mean <- round(colMeans(Beta1) ,  2)
beta1.ci <- round(apply(Beta1 ,  2 ,  quantile ,  prob = c(0.025 , 0.975)) ,  2)

####################
# count component  #
####################

beta2.mean <- round(colMeans(Beta2) ,  2)
beta2.ci <- round(apply(Beta2 ,  2 ,  quantile ,  prob = c(0.025 , 0.975)) ,  2)


#################################
# Overdispersion parameter: R   #
#################################

r.est <- round(mean(R) ,  2)
r.ci <- round(quantile(R ,  probs = c(0.025 ,  0.975)) ,  2)


#######################
# covariance matrix   #
#######################

sigmaphi.est <- round(apply(Sigphi ,  2 ,  mean) ,  2)
sigmaphi.ci <- round(apply(Sigphi ,  2 ,  quantile ,  probs = c(0.025 ,  0.975)) ,  2)



####################
# Table 1          #
####################

###################################################################################################
# A data frame for Table 1. The Greek symbols for parameters are presented in character form here #
###################################################################################################


df_table1 <- tibble(Parameter = c("gamma" ,  "beta" ,  "Sigma11" ,  "Sigma12" ,  "Sigma22" ,  "r") , 
                    Truth = c(beta1.true[1] ,  beta2.true[1] ,  cov[1 ,  1] ,  cov[1 ,  2] ,  cov[2 ,  2] ,  1) ,  
                    `Post. Mean` = c(beta1.mean[1] ,  beta2.mean[1] ,  sigmaphi.est[1] ,  sigmaphi.est[2] ,  sigmaphi.est[4] ,  r.est) , 
                    `95% LCrI` = c(beta1.ci[1 ,  1] ,  beta2.ci[1 ,  1] ,  sigmaphi.ci[1 ,  1] ,  sigmaphi.ci[1 ,  2] ,  sigmaphi.ci[1 ,  4] ,  r.ci[1]) , 
                    `95% UCrI` = c(beta1.ci[2 ,  1] ,  beta2.ci[2 ,  1] ,  sigmaphi.ci[2 ,  1] ,  sigmaphi.ci[2 ,  2] ,  sigmaphi.ci[2 ,  4] ,  r.ci[2]))

table1 <- setDT(df_table1)
table1

write_csv(table1 ,  file = "./bimj_files/results/Table_1.csv")







###############################################
# Time trends for the binary and count parts  #
###############################################

###################
#  binary part    #
###################

ft.bin.true <- c(splines.tmp%*%alpha1)
est.alpha1 <- colMeans(Beta1[ ,  -c(1)])
ft.est.bin <- c(splines.tmp%*%est.alpha1)

lcl.bin <- apply(Beta1[ ,  -c(1)] ,  2 ,  quantile ,  prob = c(0.025))
lbin_trend <- c(splines.tmp%*%lcl.bin)

ucl.bin <- apply(Beta1[ ,  -c(1)] ,  2 ,  quantile ,  prob = c(0.975))
ubin_trend <- c(splines.tmp%*%ucl.bin)


###############################
# count (marginal mean) part  #
###############################

ft.mm.true <- c(splines.tmp%*%alpha2)
est.alpha2 <- colMeans(Beta2[ ,  -c(1)])
ft.est.mm <- c(splines.tmp%*%est.alpha2)

lcl.mm <- apply(Beta2[ ,  -c(1)] ,  2 ,  quantile ,  prob = c(0.025))
lmm_trend <- c(splines.tmp%*%lcl.mm)

ucl.mm <- apply(Beta2[ ,  -c(1)] ,  2 ,  quantile ,  prob = c(0.975))
umm_trend <- c(splines.tmp%*%ucl.mm)

########################################
# A data frame for the time trends     #
########################################

splines_data <- tibble(day = 1:lday , ft.bin.true = ft.bin.true ,  ft.est.bin = ft.est.bin ,  ft.mm.true = ft.mm.true ,  ft.est.mm = ft.est.mm ,  lbin_trend = lbin_trend , 
                       ubin_trend = ubin_trend ,  lmm_trend = lmm_trend ,  umm_trend = umm_trend)


#######################################################
# Time trends for the binary component: Figure 3 (a)  #
#######################################################

splines_binary <- ggplot(splines_data ,  aes(x = day ,  y = ft.bin.true ))+
  geom_line(aes(color = 'Simulated Curve') ,  size = 0.5 ,  linetype = 1)+
  geom_line(aes(day ,  ft.est.bin ,  color = "Posterior Mean Curve") ,  size = 0.5 ,  linetype = 2)+
  geom_ribbon(aes(ymin = lbin_trend ,  ymax = ubin_trend ,  color="95% Credible Interval") , alpha = 0.1 ,  show.legend = F)+
  scale_color_manual(breaks = c("Simulated Curve" ,  
                                "Posterior Mean Curve" , 
                                "95% Credible Interval"
  ) ,  
  values = c("Simulated Curve" = "darkorange" , 
             "Posterior Mean Curve" = "blue4" , 
             "95% Credible Interval" = "gray40"
  ))+
  ylab(expression(f[1](t)))+
  ylim(-2.75 ,  2.75)+
  coord_fixed(ratio = 150/(5.5))+
  scale_x_continuous(name = "Day" , breaks = c(1 ,  32 ,  60 ,  91 ,  121 ,  152 ,  182 ,  213 ,  244 ,  274 ,  305 ,  335 ,  365)
                     #labels=c("Jan 1" ,  "Feb 1" ,  "March 1" ,  "Apr 1" ,  "May 1" ,  "June 1" , "July 1" ,  "Aug 1" ,  "Sep 1" ,  "Oct 1" ,  "Nov 1" ,  "Dec 1" ,  "Dec 31")
  )+
  guides( color = guide_legend(title = "Time Trend: Binary part" , 
                               override.aes = list(
                                 linetype = c(1 ,  2 ,  1) , 
                                 shape = c(16 ,  NA ,  NA)) , 
                               reverse = F))+
  ggtitle("(a)")+
  theme(legend.key = element_rect(fill = "white") , 
        legend.position = c(0.4 ,  0.8) , 
        legend.text = element_text(size = 8) , 
        legend.title = element_text(size = 8 ,  face = 'bold') , 
        axis.text = element_text(size = 8 ,  face = "bold") , 
        axis.title = element_text(size = 8) , 
        plot.title = element_text(size = 8 ,  hjust = 0.5 ,  face = "bold") , 
        panel.background = element_rect(colour = "grey50") , 
        plot.margin = unit(c(0 , 0 , 0 , 0) ,  "cm"))


ggsave(
  "./bimj_files/results/Figure_3a.pdf" , 
  plot = splines_binary , 
  device = pdf() , 
  scale = 1 , 
  dpi = "print"(300) , 
  width = 7 , 
  height = 3.5 , 
  units = "in"
)

dev.off()


##################################################
# Time trends for the count part: Figure 3 (b)   #
##################################################

splines_count <- ggplot(splines_data ,  aes(x = day ,  y = ft.mm.true ))+
  geom_line(aes(color = 'Simulated Curve') ,  size = 0.5 ,  linetype = 1)+
  geom_line(aes(day ,  ft.est.mm ,  color = "Posterior Mean Curve") ,  size = 0.5 ,  linetype = 2)+
  geom_ribbon(aes(ymin = lmm_trend ,  ymax = umm_trend ,  color = "95% Credible Interval") , alpha = 0.1 ,  show.legend = F)+
  scale_color_manual(breaks = c("Simulated Curve" ,  
                                "Posterior Mean Curve" , 
                                "95% Credible Interval"
  ) ,  
  values = c("Simulated Curve" = "darkorange" , 
             "Posterior Mean Curve" = "blue4" , 
             "95% Credible Interval" = "gray40"
  ))+
  ylab(expression(f[2](t)))+
  ylim(-2.25 ,  1.5)+
  coord_fixed(ratio = 150/(3.75))+
  scale_x_continuous(name="Date" , breaks = c(1 ,  32 ,  60 ,  91 ,  121 ,  152 ,  182 ,  213 ,  244 ,  274 ,  305 ,  335 ,  365)
                     #labels = c("Jan 1" ,  "Feb 1" ,  "March 1" ,  "Apr 1" ,  "May 1" ,  "June 1" , "July 1" ,  "Aug 1" ,  "Sep 1" ,  "Oct 1" ,  "Nov 1" ,  "Dec 1" ,  "Dec 31")
  )+
  guides( color = guide_legend(title = "Time Trend: Overall mean" , 
                               override.aes = list(
                                 linetype = c(1 ,  2 ,  1) , 
                                 shape = c(NA ,  NA ,  NA)) , 
                               reverse = F))+
  ggtitle("(b)")+
  theme(legend.key = element_rect(fill = "white") , 
        legend.position = c(0.4 ,  0.8) , 
        legend.text = element_text(size = 8) , 
        legend.title = element_text(size = 8 ,  face = 'bold') , 
        axis.text = element_text(size = 8 ,  face = "bold") , 
        axis.title = element_text(size = 8) , 
        plot.title = element_text(size = 8 ,  hjust = 0.5 ,  face = "bold") , 
        panel.background = element_rect(colour = "grey50") , 
        plot.margin = unit(c(0 , 0 , 0 , 0) ,  "cm"))

ggsave(
  "./bimj_files/results/Figure_3b.pdf" , 
  plot = splines_count , 
  device = pdf() , 
  scale = 1 , 
  dpi = "print"(300) , 
  width = 7 , 
  height = 3.5 , 
  units = "in"
)

dev.off()



#############################
# Random effect maps        #
#############################

ga_counties_map <- counties(state = c('13'))

GEOID <- ga_counties_map$GEOID



############################################
# Simulated vs Predicted estimate data     #
############################################

mean_phi1 <- apply(PHI1 ,  2 ,  mean)
mean_phi2 <- apply(PHI2 ,  2 ,  mean)

re_means <- tibble(true_phi1 = c(true.phi1) ,  mean_phi1 = c(mean_phi1) ,  true_phi2 = c(true.phi2) ,  mean_phi2 = c(mean_phi2) ,  GEOID = GEOID)
re_data <- left_join(ga_counties_map ,  re_means ,  by = "GEOID")   # merge the random effects estimates data with county level geoid data


pal <- brewer.pal(5 , "BuGn")        # specify palette colors


#######################################################
# map for the binary component: Figures 4 (a) and (b) #
#######################################################

tm_binary <- tm_shape(re_data)+
  tm_fill(c("true_phi1" ,  "mean_phi1") ,  midpoint = c(NA) ,  title = c(expression(paste("Simulated "~phi[1])) ,  expression(paste("Predicted"~phi[1]))) ,  palette = pal ,  style = "quantile")+
  tm_facets(free.scales = FALSE ,  nrow = 1)+
  tm_layout(title = "Quintile Map" , 
            title.snap.to.legend = TRUE , 
            panel.labels = c("(a)" ,  "(b)") , 
            panel.label.bg.color = c('white') , 
            panel.label.fontface = c('bold') , 
            title.size = 0.8 , 
            title.position = c("right" ,  "bottom") , 
            legend.outside = FALSE , 
            legend.position = c(0.65 ,  0.75))+
  tm_borders(alpha = 0.3 ,  lwd = 1)


tmap_save(tm_binary ,  filename = "./bimj_files/results/Figures_4_ab.jpg" ,  dpi = 300 )


######################################################################
# map for the count (marginal mean) component: Figures 4 (c) and (d) #
######################################################################

tm_count <- tm_shape(re_data)+
  tm_fill(c("true_phi2" ,  "mean_phi2") ,  midpoint = c(NA) ,  title = c(expression(paste("Simulated "~phi[2])) ,  expression(paste("Predicted"~phi[2]))) ,  palette = pal ,  style = "quantile")+
  tm_facets(free.scales = FALSE ,  nrow = 1)+
  tm_layout(title = "Quintile Map" , 
            title.snap.to.legend = TRUE , 
            panel.labels = c("(c)" ,  "(d)") , 
            panel.label.bg.color = c('white') , 
            panel.label.fontface = c('bold') , 
            title.size = 0.8 , 
            title.position = c("right" ,  "bottom") , 
            legend.outside = FALSE , 
            legend.position = c(0.65 ,  0.75))+
  tm_borders(alpha = 0.3 ,  lwd = 1)


tmap_save(tm_count ,  filename = "./bimj_files/results/Figures_4_cd.jpg" ,  dpi = 300 )

