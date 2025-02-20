##############################################################################################################################
#                                                                                     
#   Filename    :    case_study.R    												  
#                                                                                     
#   Project     :    BiomJ article "A marginalized zero-inflated negative binomial
#                    model for spatial data: modeling COVID-19 deaths in Georgia"   
#   Authors     :    Fedelis Mutiso , Hong Li , John L. Pearce , Sara E. Benjamin-Neelon , Noel T. Mueller , and Brian Neelon                                                              
#   Date        :    05.20.2024
#   Purpose     :    Reproduce Figure_1.jpg , Figure_2.pdf , Table_2.csv , Figure_5.pdf , and Figure_6.jpg of the manuscript
#																				  
#   R Version   :    4.2.2 (2022-10-31 ucrt)                                                                
#
#   Input data files  :    covid_data_ga.csv , ordered_adjmat_ga.txt , and beta.inits.csv  
#
#   Output data files :    Figure_1.jpg , Figure_2.pdf , Table_2.csv , Figure_5.pdf , Figure_6.jpg
#
################################################################################################################################



library(BayesLogit)   # For rpg function 
library(mvtnorm)	
library(MCMCpack)     # For Iwish update of Sigmab
library(msm)          # For tnorm function
library(spam)         # For sparse Matrix
library(splines2)
library(tidyverse)
library(pscl)
library(data.table)

library(tigris)
library(rgdal)
library(sf)
library(tmap)
library(RColorBrewer)


setwd('~/')                 # specify the full path to the location where the `bimj_files` folder containing the r scripts and data/results subfolders is stored


set.seed(091322)


##############################
# read in the data           #
##############################

covid_data_ga <- read.csv(file = './bimj_files/data/covid_data_ga.csv' , header = TRUE)

fips <- covid_data_ga[ , 1]
ncounty <- length(unique(fips))
lday <- nrow(covid_data_ga)/ncounty
days <- rep(1:lday , ncounty)
nis <- rep(lday , ncounty)

attach(covid_data_ga)


#################################
# Adjacency matrix information  #
#################################

A <- matrix(scan("./bimj_files/data/ordered_adjmat_ga.txt") , 159 , 159)
m <- apply(A , 2 , sum)
kappa <- .999999			  	                      # Spatial dependence parameter ~ 1 for intrinsic CAR
Q <- as.spam(diag(m))-kappa*as.spam(A)      

adjid <- rep(1:ncounty , m)                    # Tells which elements of A belong to which block
adj <- rep(0 , sum(m))                         # how many in total are adjacent
for (i in 1:ncounty) adj[adjid==i] <- which(A[i , ]==1)


######################################################################
# Adjustment Variables from RWJ Foundation , US CENSUS , NOAA , or CDC  # 
######################################################################

svi <- scale(svi_overall_rank)
smoke <- scale(pct_adult_smoke)
phys <- scale(no_of_pc_phys_per100k)   
poor_health <- scale(pct_poor_or_fair_health)
cvd_hosp <- scale(cvd_hosp_rate)
pop_dens <- scale(popn_density)
pm25 <- scale(average_daily_PM_25)
temp <- scale(average_temp)
precip <- scale(precipitation)



###################################################################################################################
# Impute any missing values with the mean of the non-missing values. However , data for most variables is complete #
###################################################################################################################


svi[is.na(svi)==T] <- mean(svi , na.rm = T)
smoke[is.na(smoke)==T] <- mean(smoke , na.rm = T)
phys[is.na(phys)==T] <- mean(phys , na.rm = T) 
poor_health[is.na(poor_health)==T] <- mean(poor_health , na.rm = T)
cvd_hosp[is.na(cvd_hosp)==T] <- mean(cvd_hosp , na.rm = T)
pop_dens[is.na(pop_dens)==T] <- mean(pop_dens , na.rm = T)
pm25[is.na(pm25)==T] <- mean(pm25 , na.rm = T)
temp[is.na(temp)==T] <- mean(temp , na.rm = T)
precip[is.na(precip)==T] <- mean(precip , na.rm = T)



lpop <- log(Population)                                                   # county population sizes


##################################################################################
# Figure 1: Descriptive map for cumulative incidence death rates                 #
##################################################################################


###################################################
# create a data frame for total deaths by county  #
###################################################

total_deaths <- covid_data_ga%>%
                select(FIPS , Incident_Deaths , Population)%>%
                mutate(GEOID = as.character(FIPS))%>%
                group_by(GEOID , Population)%>%
                summarise(n_deaths = sum(Incident_Deaths))%>%
                mutate(cum_ir = n_deaths/Population*1e6)



ga_counties_map <- counties(state = c('13'))

deaths_data <- left_join(ga_counties_map , total_deaths , by = 'GEOID')

pal <- brewer.pal(5 , "BuGn")  # specify the palette colors


########################################################################
# Figure 1 : A map of crude cumulative COVI-19 death rates per million #
########################################################################


figure1 <- tm_shape(deaths_data)+
            tm_fill(c("cum_ir") , midpoint = c(NA) , title = c(expression(paste("Deaths per Million"))) , palette = pal , style = "quantile")+
            tm_facets(free.scales = TRUE , nrow = 1)+
            tm_layout(title = "Quintile Map" , 
                      title.snap.to.legend = TRUE , 
                      title.size = 0.8 , 
                      title.position = c("right" , "bottom") , 
                      legend.outside = FALSE , 
                      legend.position = c(0.7 , 0.85) , 
                      main.title.fontface = "bold" , 
                      main.title.position = "center")+
            tm_borders(alpha = 0.3 , lwd = 1)


tmap_save(figure1 , filename = "./bimj_files/results/Figure_1.jpg" , dpi = 300 )


############################################
# Figure 2 : A bar plot of COVID-19 deaths #
############################################

histdata <- tibble(y = Incident_Deaths)

bar_plot <- ggplot(histdata)+
            geom_bar(aes( x = y , y = ..prop.. , group = 1) , stat = "count")+
            scale_y_continuous(labels = scales::label_percent(accuracy = 1L))+
            ylab("Percent")+
            xlab("Number of Deaths")


ggsave(
  "./bimj_files/results/Figure_2.pdf" , 
  plot = bar_plot , 
  device = pdf() , 
  scale = 1 , 
  dpi = "print"(300) , 
  width = 6 , 
  height = 6 , 
  units = "in"
)

dev.off()




##################
# Bayesian Model #
##################


knots <- seq(14 , 351 , by = 14)  
covariates <- cbind(1 , svi , smoke , phys , poor_health , cvd_hosp , pop_dens , pm25 , temp , precip)                                                    
Xtmp <- bSpline(1:lday , knots = knots , intercept = F)
head(Xtmp)
k <- ncol(Xtmp)                                           # Number of spine coefs
U <- Xspline <- apply(Xtmp , 2 , rep , ncounty)
X1 <- X2 <- X <- cbind(covariates , Xspline)
n <- nrow(X)
p <- ncol(X)
nis <- rep(lday , ncounty)
id <- rep(1:ncounty , nis)
N <- length(id)



##########   
# Priors #
##########

beta.inits <- read_csv(file = './bimj_files/data/beta.inits.csv' , col_names = TRUE)   

beta10 <- c(beta.inits$beta1)
beta20 <- c(beta.inits$beta2)
# beta10 <- rep(0 , p)
# beta20 <- rep(0 , p)
T0b1 <- diag(.01 , p)
s <- 0.015



####################################
# Initialize different parameters  #
####################################

# beta1 <- rep(0 , p)
# beta2 <- rep(0 , p)

beta1 <- beta10
beta2 <- beta20

phi_init <- spam::rmvnorm(1 , sigma = diag(.1 , 2*ncounty))	  # Random effects
phi_init <- matrix(phi_init , ncol = 2 , byrow = T)               # n x 2 matrix of spatial effects
phi1 <- phi_init[ , 1]
phi2 <- phi_init[ , 2]
Phi1 <- rep(phi1 , nis)
Phi2 <- rep(phi2 , nis)


phibar1 <- phibar2 <- rep(0 , ncounty)      # mean of adjacent phi2's for updating phi2
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

bcov <- diag(0.1 , p)
sb2 <- 0.005               # tuning parameter for covariance of parameters for the mean component
accb <- 0


s2 <- 0.17                  # tuning parameter for the variance of the random effects for the mean
u <- rep(0 , ncounty)       # keeps track of acceptance of phi2 for each county at each iteration
accu <- 0                   # keeps track of the mean acceptance rate of phi2

y <- Incident_Deaths
r <- 1
accr <- 0  
y1 <- rep(0 , N)             # At risk indicator (this is W in paper)
y1[y>0] <- 1                 # If y>0 , then at risk w.p. 1
N0 <- length(y[y==0])        # Number of observed 0's
q <- rep(.5 , N)             # 1-p=1/(1+exp(X%*%alpha)) , used for updating y1


############
# Num Sims #
############
nsim <- 100000			        # Number of MCMC Iterations (100k)
thin <- 25		                # Thinning interval (25)
burn <- 50000		                # Burnin (50k)
lastit <- (nsim-burn)/thin	        # Last stored value

#########
# Store #
#########
Beta1 <- matrix(0 , lastit , p)
Beta2 <- matrix(0 , lastit , p)
# BETA1 <- matrix(0 , nsim , p)             
BETA2 <- matrix(0 , nsim , p)             # for updating beta2 proposal
R <- R2 <- rep(0 , lastit)
Sigphi <- matrix(0 , lastit , 4)
PHI1 <- PHI2 <- matrix(0 , lastit , ncounty)


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
  mb <- vb%*%(T0b1%*%beta10+t(sqrt(w)*X1)%*%(sqrt(w)*(z-Phi1))) #-Phi1-Phi2*x
  beta1 <- c(spam::rmvnorm(1 , mb , vb))
  
  # BETA1[i , ] <- beta1
  
  # update phi1
  
  priorprec <- 1/(sphi1^2*(1-rho^2))*Q                                               # Prior precision of phi1|phi2
  priormean <- rho*sphi1/sphi2*(phi2)                                                # Prior mean of phi1|phi2 --- 
  prec <- priorprec+as.spam(diag(tapply(w , id , sum) , ncounty , ncounty))
  mb <- c(priorprec%*%priormean)+tapply(w*(z-X1%*%beta1) , id , sum)                         
  phi1 <- rmvnorm.canonical(1 , mb , prec)[1 , ]
  
  
  # center phi1 and update tauphi1

  phi1 <- phi1-mean(phi1)
  Phi1 <- rep(phi1 , nis)
  
  # update phibar   (apply the sum to zero constraint)

  for (j in 1:ncounty) phibar1[j] <- mean(phi1[adj[which(adjid==j)]])
  
  
  ###################################
  # update the at-risk indicator    #
  ###################################
  
  
  eta1 <- X1%*%beta1+Phi1
  eta2 <- X2%*%beta2+lpop-log(10^6)+Phi2                          # Use all n observations
  psi <- c(pmax(0.01 , pmin(0.99 , exp(eta1)/(1+exp(eta1)))))     # at-risk probability
  nu = c(exp(eta2))
  mu <- nu/psi
  q <- r/(mu+r)                                    # Pr(y=0|y1=1)
  theta <- psi*(q^r)/(psi*(q^r)+1-psi)             # Conditional prob that y1=1 given y=0 -- i.e. Pr(chance zero|observed zero)
  y1[y==0] <- rbinom(N0 , 1 , theta[y==0])         # If y=0 , then draw a "chance zero" w.p. theta , otherwise y1=1
  N1 <- sum(y1)
  nis1 <- tapply(y1 , id , sum)
  
  ###############
  # Update r    #
  ###############
  
  rnew <- rtnorm(1 , r , sqrt(s) , lower = 0)       # Treat r as continuous
  ratio <- sum(dnbinom(y[y1==1] , rnew , mu = mu[y1==1] , log = T))-sum(dnbinom(y[y1==1] , r , mu = mu[y1==1] , log = T))+   
    dtnorm(r , rnew , sqrt(s) , 0 , log = T) - dtnorm(rnew , r , sqrt(s) , 0 , log = T)       # diffuse prior for r
  # Proposal not symmetric 
  if (log(runif(1))<ratio) {
    r <- rnew
    accr <- accr+1
  }
  
  
  ################
  # update beta2 #
  ################
  
  # Proposal
  beta2s <- c(beta2+rmvt(1 , df = 3 , sigma = sb2*bcov))
  eta2s <- X2[y1==1 , ]%*%beta2s+lpop[y1==1]-log(10^6)+Phi2[y1==1]
  nus <- c(exp(eta2s))
  mus <- nus/psi[y1==1]
  ls <- sum(dnbinom(y[y1==1] , r , mu = mus , log = T))                                 # Log like under proposal
  lps <- dmvnorm(beta2s , mean = rep(0 , p) , sigma = diag(10 , p) , log = T)           # Log prior proposal
  
  # Current
  eta2 <- X2[y1==1 , ]%*%beta2+lpop[y1==1]-log(10^6)+Phi2[y1==1]
  nu <- c(exp(eta2))
  mu <- nu/psi[y1==1]
  l <- sum(dnbinom(y[y1==1] , r , mu = mu , log = T))
  lp <- dmvnorm(beta2 , mean = rep(0 , p) , sigma = diag(10 , p) , log = T)
  
  ratio <-  ls+lps-l-lp    
  
  if (log(runif(1))<ratio) {
    beta2 <- beta2s
    accb <- accb+1
  }
  
  BETA2[i , ] <- beta2
  if (i==nsim/2) bcov <- cov(BETA2[(nsim/4+1):(nsim/2) , ])  # Tune proposal covariance for Beta
  
  
  
  ###############
  # Update phi2 #
  ###############
  
  # likelihood - Proposal
  
  phi2s <- phi2+s2*rt(ncounty , 3)
  Phi2s <- rep(phi2s , nis)
  eta2s <- X2[y1==1 , ]%*%beta2+Phi2s[y1==1]+lpop[y1==1]-log(10^6)
  nus <- c(exp(eta2s))
  mus <- nus/psi[y1==1]
  Ls <- dnbinom(y[y1==1] , r , mu = mus , log = T)
  ls <- rep(0 , ncounty)                              # account for empty counties - their likelihood contribution is 0
  ls[nis1>0] <- tapply(Ls , id[y1==1] , sum)
  
  # likelihood - Current
  eta2 <- X2[y1==1 , ]%*%beta2+Phi2[y1==1]+lpop[y1==1]-log(10^6)
  nu <- c(exp(eta2))
  mu <- nu/psi[y1==1]
  L <- dnbinom(y[y1==1] , r , mu = mu , log = T)
  l <- rep(0 , ncounty)
  l[nis1>0] <- tapply(L , id[y1==1] , sum)
  
  # prior -- both current and proposal
  for (j in 1:ncounty){ 
    lp <- dnorm(phi2[j] , phibar2[j]+rho*sphi2/sphi1*(phi1[j]-phibar1[j]) , sqrt((sphi2^2*(1-rho^2)/m[j])) , log = T)
    lps <- dnorm(phi2s[j] , phibar2[j]+rho*sphi2/sphi1*(phi1[j]-phibar1[j]) , sqrt((sphi2^2*(1-rho^2)/m[j])) , log = T)
    ratio <- ls[j]+lps-l[j]-lp
    u[j] <- 1*(log(runif(1))<ratio)                    # keep track of acceptance of phi2 for each county
    if (log(runif(1))<ratio) {
      phi2[j] <- phi2s[j]
    }
    phibar2[j] <- mean(phi2[adj[which(adjid==j)]])     # only updated if current phi2 is accepted
  }
  
  accu[i] <- mean(u)
  
  # Center phi2
  phi2 <- phi2-mean(phi2)
  Phi2 <- rep(phi2 , nis)
  
  # Update Sigma.phi
  
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
    PHI1[j , ] <- c(phi1)
    PHI2[j , ] <- c(phi2)
    
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




#############################
# posterior mean estimates  #
#############################


# Binary component

beta1_fx_means <- round(colMeans(Beta1)[1:10] , 2)                           # only the fixed effects
beta1_fx_ci <- round(apply(Beta1 , 2 , quantile , prob = c(0.025 , 0.975))[ , 1:10] , 2)


# Count part

beta2_fx_means <- round(colMeans(Beta2)[1:10] , 2)                          # only the fixed effects
beta2_fx_ci <- round(apply(Beta2 , 2 , quantile , prob = c(0.025 , 0.975))[ , 1:10] , 2)


# NB overdispersion parameter r

r.mean <- round(mean(R) , 2)
r.ci <- round(quantile(R , probs = c(0.025 , 0.975)) , 2)

# covariance matrix

Sigmaphi.est <- round(apply(Sigphi , 2 , mean) , 2)
Sigmaphi.ci <- round(apply(Sigphi , 2 , quantile , probs = c(0.025 , 0.975)) , 2)



########################
# Table 2              #
########################

###################################################################################################
# A data frame for Table 2. The Greek symbols for parameters are presented in character form here #
###################################################################################################

var_bin <- c("SVI" , "% of adult smokers" , "No. of physicians per 100K" , "% fair or poor health" , 
             "CVD Hospitalizations" , "Population density" , "PM2.5" , "Temperature" , "Precipitation")
var_re <- c("var_phi1i" , "cov_phi1i_phi2i" , "var_phi2i")

df_table2 <- tibble(`Model Component` = c("Binary" , rep('' , 9) , "Mean" , rep('' , 10) , "Random Effects" , rep('' , 2)) , 
                    Variable = c(var_bin , '' , var_bin , 'Dispersion' , '' , var_re) , 
                    Parm = c("gamma1" , "gamma2" , "gamma3" , "gamma4" , "gamma5" , "gamma6" , "gamma7" , "gamma8" , "gamma9" , "" , 
                             "beta1" , "beta2" , "beta3" , "beta4" , "beta5" , "beta6" , "beta7" , "beta8" , "beta9" , "r" , "" , 
                             "Sigma11" , "Sigma12" , "Sigma22") , 
                    `Posterior Mean` = c(beta1_fx_means[-1] , NA , beta2_fx_means[-1] , r.mean , NA , Sigmaphi.est[-3]) , 
                    `95% LCrI` = c(beta1_fx_ci[1 , -1] , NA , beta2_fx_ci[1 , -1] , r.ci[1] , NA , Sigmaphi.ci[1 , -3]) , 
                    `95% UCrI` = c(beta1_fx_ci[2 , -1] , NA , beta2_fx_ci[2 , -1] , r.ci[2] , NA , Sigmaphi.ci[2 , -3]))

table2 <- setDT(df_table2)
table2

write_csv(table2 , file = "./bimj_files/results/Table_2_R1.csv")





######################
# Figure 5
######################


#############################################
# Time trends estimates
############################################

Pred <- array(0 , dim=c(ncounty , lday , lastit))
for (j in 1:lastit){
  for (h in 1:ncounty){
    beta <- Beta2[j , 11:38]
    alpha <- Beta2[j , 1:10]
    Pred[h , , j] <- c(exp(Xtmp%*%beta+alpha[1]*mean(covariates[id==h , 1])+alpha[2]*mean(covariates[id==h , 2])
                      +alpha[3]*mean(covariates[id==h , 3])+alpha[4]*mean(covariates[id==h , 4])+alpha[5]*mean(covariates[id==h , 5])
                      +alpha[6]*mean(covariates[id==h , 6])+alpha[7]*mean(covariates[id==h , 7])+alpha[8]*mean(covariates[id==h , 8])
                      +alpha[9]*mean(covariates[id==h , 9])+alpha[10]*mean(covariates[id==h , 10])))
  } 
  if (j%%100==0){
    print(j)
  }
}


ft_mzinb <- predMean <- colMeans(t(colMeans(Pred[ , , 1:lastit])))     # population average prediction (average across counties and iterations)
lcl_mzinb <- palower_mean <- apply(t(colMeans(Pred[ , , 1:lastit])) , 2 , quantile , prob = 0.025)
ucl_mzinb <- paupper_mean <- apply(t(colMeans(Pred[ , , 1:lastit])) , 2 , quantile , prob = 0.975)

day <- c(1:lday)

incidence <- y/Population*1e6
orr <- tapply(incidence , days , mean)                              # average daily incidence across all counties
splines_data <- tibble(day = 1:lday , orr = orr , ft.est = predMean , lcl = palower_mean , ucl = paupper_mean)
splines_mzinb <- as.data.frame(cbind(day , orr , ft_mzinb , lcl_mzinb , ucl_mzinb))


figure5 <- ggplot(splines_data , aes(x = day , y = orr ))+
            geom_point(aes(color = 'Daily Average') , shape = 16 , size = 1)+
            geom_line(aes(day , ft.est , color = "Estimated Trend") , size = 0.5 , linetype = 1)+
            geom_ribbon(aes(ymin = lcl , ymax = ucl , color="95% Credible Interval") , alpha = 0.1 , show.legend = F)+
            scale_color_manual(breaks = c("Daily Average" , 
                                          "Estimated Trend" , 
                                          "95% Credible Interval"
            ) , 
            values = c("Daily Average" = "darkorange" , 
                       "Estimated Trend" = "blue4" , 
                       "95% Credible Interval" = "gray40"
            ))+
            ylab('Deaths per million')+
            ylim(0 , 40)+
            coord_fixed(ratio = 150/(40))+
            scale_x_continuous(name = "Date" , breaks = c(1 , 32 , 60 , 91 , 121 , 152 , 182 , 213 , 244 , 274 , 305 , 335 , 365) , 
                               labels = c("Jan 1" , "Feb 1" , "March 1" , "Apr 1" , "May 1" , "June 1" , "July 1" , "Aug 1" , "Sep 1" , "Oct 1" , "Nov 1" , "Dec 1" , "Dec 31")
            )+
            guides( color = guide_legend(title = "Fatality Rate Trends" , 
                                         override.aes = list(
                                           linetype = c(0 , 1 , 1) , 
                                           shape = c(16 , NA , NA)) , 
                                         reverse = F))+
            # ggtitle("()")+
            theme(legend.key = element_rect(fill = "white") , 
                  legend.position = c(0.4 , 0.8) , 
                  legend.text = element_text(size = 8) , 
                  legend.title = element_text(size = 8 , face = 'bold') , 
                  axis.text = element_text(size = 8 , face = "bold") , 
                  axis.title = element_text(size = 8) , 
                  plot.title = element_text(size = 8 , hjust = 0.5 , face = "bold") , 
                  panel.background = element_rect(colour = "grey50") , 
                  plot.margin = unit(c(0 , 0 , 0 , 0) , "cm"))



ggsave(
  "./bimj_files/results/Figure_5_R1.pdf" , 
  device = pdf() , 
  plot = figure5 , 
  scale = 1 , 
  dpi = "print"(300) , 
  width = 7 , 
  height = 3.5 , 
  units = "in"
)

dev.off()




#########################################################
# Figure 6: Random effect maps for the count component  #
#########################################################


ga_counties_map <- counties(state = c('13'))

GEOID <- ga_counties_map$GEOID


mean_phi2 <- apply(PHI2 , 2 , mean)

re_means <- tibble(mean_phi2 = c(mean_phi2) , GEOID = GEOID)
re_data <- left_join(ga_counties_map , re_means , by = "GEOID")   # merge the random effects estimates data with county level geoid data


pal <- brewer.pal(5 , "BuGn")        # specify palette colors

figure6 <- tm_shape(re_data)+
            tm_fill(c("mean_phi2") , midpoint = c(NA) , title = c(expression(paste("Predicted"~phi[2]))) , palette = pal , style = "quantile")+
            tm_facets(free.scales = FALSE , nrow = 1)+
            tm_layout(title = "Quintile Map" , 
                      title.snap.to.legend = TRUE , 
                      panel.label.bg.color = c('white') , 
                      panel.label.fontface = c('bold') , 
                      title.size = 0.8 , 
                      title.position = c("right" , "bottom") , 
                      legend.outside = FALSE , 
                      legend.position = c(0.65 , 0.75))+
            tm_borders(alpha = 0.3 , lwd = 1)


tmap_save(figure6 , filename = "./bimj_files/results/Figure_6_R1.jpg" , dpi = 300 )

