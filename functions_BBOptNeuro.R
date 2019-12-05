
#Functions for Batch Bayesian Optimization Design for Optimizing a Neurostimulator
# By Adam Kaplan and Thomas Murray (2019)
# University of Minnesota, Twin Cities 

#Download these prior #
#install.packages("emdbook")
#install.packages("INLA")
#install.packages("BayesLogit")

require(emdbook)
require(INLA)
library(BayesLogit)
############################################################################

###Function to generate a random ordering of the 8 parameter values within each month
#Ensures the same parameter setting is not repeated on consequtive days
#Ensures each setting is given the maximum possible times, e.g., 3 or 4 times for 8 settings over 30 days
get.order = function(nvalues,ndays){
  #Generate random ordering
  order = sample(rep(1:nvalues,ceiling(ndays/nvalues)),ndays)
  
  #Re-sample the ordering until no setting is repeated on consequtive days
  while(sum(order[-ndays]==order[-1])) order = sample(rep(1:nvalues,ceiling(ndays/nvalues)),ndays)
  
  return(order)
}

# Simulates 1-D True Latent Surfaces #
estBetaParams <- function(mu, var) {
  alpha <- mu*((mu*(1-mu))/var - 1)
  beta <- (1-mu)*(alpha/mu)
  return(params = list(alpha = alpha, beta = beta))
}


# another way to simulate 1-D preferences for configurations
#prevals = parameters settings (ex. 10:1000 by tens)
#cores = how many modes you want for the true latent preference curve
#mus = a vector of where the modes you desire (should be between (0, 1))
#vars = width for each hill, I highly recommend multiple hills should have 0 < vars < 0.18
#vars (continued): for one hill the max variance can be is 0.15
#and p = the mixing probability, default is 1 for unimodal 
pref.funcval <- function(prefvals, cores, mus, vars, p = 1, Qtype = 1){
  input <- prefvals/max(prefvals); 
  out <- numeric(length(prefvals)) 
  for(i in 1:cores){
    conv <- estBetaParams(mus[i], vars[i])
    a = conv$alpha; b = conv$beta
    
    output <- p[i]*dbeta(input, a, b);
    out <- out + output
  }
  #need to fix this when things are infinite on the tails #
  out <- out[is.finite(out)]
  
  output2 <- scale(out) #does this center the data?
  output2 <- as.numeric(output2);
  #plot(output2~prefvals, xlab = "Parameter Values", ylab = "Latent Preference", main = "True Latent Surface", type = "l")
  return(output2)
}

estBetaParams2 <- function(mu, sumab) {
  alpha <- mu * (sumab - 2) + 1
  beta <- sumab - alpha
  return(params = list(alpha = alpha, beta = beta))
}

pref.funcvalSUM <- function(prefvals, cores, mus, sumab, p = 1){
  input <- prefvals/max(prefvals); 
  out <- numeric(length(prefvals)) 
  for(i in 1:cores){
    conv <- estBetaParams2(mus[i], sumab[i])
    a = conv$alpha; b = conv$beta
    
    output <- p[i]*dbeta(input, a, b);
    out <- out + output
  }
  #need to fix this when things are infinite on the tails #
  out <- out[is.finite(out)]
  
  output2 <- scale(out) #does this center the data?
  output2 <- as.numeric(output2);
    return(output2)
}

###Function to generate one month's data (i.e. indicators for preference between current and previous setting)
get.month.data2 = function(alphatrue, order, batch){
  #Preference probability each day is the expit of the difference between latent preference for the settings on current versus previous day
  
  pref.probs = 1/(1+exp(-(alphatrue[batch[order[-1]]] - alphatrue[batch[order[c(-length(order))]]])))
  #Sample preference indicators based on the calculated preference probabilities
  
  y = rbinom(length(pref.probs), 1, pref.probs)
  #Create and return dataset for the month
  
  month.data = data.frame(day=2:length(order), y=y, curr = batch[order[-1]], prev = batch[order[-length(order)]])
  return(month.data)
}

###Function to generate one month's data for 2-D configuration (i.e. indicators for preference between current and previous setting)

get.month.data2D = function(alphatrue, order, batch){
  #Preference probability each day is the expit of the difference between latent preference for the settings on current versus previous day
  
  pref.probs = 1/(1+exp(-(alphatrue[batch[order[-1]]] - alphatrue[batch[order[c(-length(order))]]])))
  y = rbinom(length(pref.probs), 1, pref.probs)
  #Sample preference indicators based on the calculated preference probabilities
  
  month.data = data.frame(day=2:length(order), y=y, curr = batch[order[-1]], prev = batch[order[-length(order)]])
  #Create and return dataset for the month
  
  return(month.data)
}

# Function that creates a rectangular lattice with indices for 2-D configuration space 
referenceLattice <- function(p1, p2){
  reffer <- matrix(1:(p1*p2), ncol = p2, nrow = p1, byrow = T)
  return(reffer)
}

# Function that converts 2D configurations to 1D monthly schedules and patient reported outcomes
# monthdat: monthly data with 1 column for outcome and 2 columns for two dimensional configuration 
# of the previous day, and 2 more columns for two dimensional configuration of the current day
# reflat: output from referenceLattice 
Convert2Dto1D <- function(monthdat, reflat){
  
  #Current Data
  matchups <- cbind(monthdat$x1.curr, monthdat$x2.curr, monthdat$x1.prev, monthdat$x2.prev)
  
  # Uses the referenceLattice function's output named "reflat" to convert D by 4 matrix
  # to a D days of matchups by 2 column matrix 
  OneDDat <- t(sapply(1:nrow(monthdat), function(x)cbind(reflat[matchups[x,3],matchups[x,4]], reflat[matchups[x,1],matchups[x,2]])))
  
  month1D <- as.data.frame(cbind(monthdat$y, OneDDat))
  names(month1D) <- c("y", "prev", "curr")
  return(month1D);
}


################################################################
# If Qtype = 2 is specified in the btgmrf function call, then 
# metropolis hasting's is used to sample for the parameter phi
# the following four functions are needed for sampling phi

# If Qtype = 1, then phi = .99 (as specified in the paper)

# Defining Q matrix for rejection sampler 
pluginQ = function(p, phi){
  Q = matrix(0,p,p); Q[1,1] = Q[p,p] = 1
  for(t in 2:(p-1)){ Q[t,t] = 1+phi^2; Q[t,t+1] = Q[t,t-1] = -phi}; Q[1,2] <- -phi; Q[p, p-1]<- -phi
  return(Q)}

#Creates iCAR Precision Matrix for one dimensional#
iCARQ = function(p){
  Q = matrix(0,p,p); counts <- numeric(p); nums <- 1:p; 
  for(t in 2:(p-1)){Q[t,t+1] = Q[t,t-1] = -1}; Q[1,2] <- Q[2,1] <- Q[p-1,p] <- Q[p, p-1] <- -1;
  counts <- sapply(nums, function(x) ((x-1) > 0) + ((x + 1) > 0)); counts[p] <- 1;
  diag(Q) = counts
  return(Q)}

#Log Posterior for Q given Beta(1,1) prior for phi#
logPostPhi2 = function(phi, alpha, tau){
  p = length(alpha); 
  val = .5*log(1-phi^2) - (p/2)*log(tau) - t(alpha)%*%pluginQ(p, phi)%*%alpha*.5; 
  return(val)
}

#Have a rejection sampler with a random walk
#The random walk adds a random uniform(0,1) to the previous value for phi#
Qreject = function(phi, alpha, tau){
  propPhi = phi + runif(1,-1,1)
  while (abs(propPhi) >= 1){propPhi <- phi + runif(1,-1,1)} #safe guarding against problematic values of phi
  if(logPostPhi2(phi = propPhi, alpha = alpha, tau = tau) - logPostPhi2(phi = phi, alpha = alpha, tau = tau) > log(runif(1))){
    newPhi <- propPhi}else{newPhi <- phi}
  return(newPhi)
}

#################################################

###Function to sample from the posterior of the Bradley-Terry GMRF Model

# par.values: configuration indices, if there are 100 configurations, then 
# par.values equals a vector of 1,...,100 (i.e., c(1:100))
# curr.data: the accumulated monthly data. One column are the binary outcomes Y's for the patient-
# reported preferences between the consecutively tested configurations, provided in the next two 
# columns "curr" and "prev". 
# an example of curr.data will be provided 

btgmrf = function(par.values, curr.data,nsamps=4000,warmup=200, Qtype = 1, hyp = c(5/2, 5/2), lambda = 0.0001){
  n = nrow(curr.data); p = length(par.values); y = curr.data$y
  a = hyp[1]; b = hyp[2]
  # Design matrix (n x p)
  X = matrix(0,nrow=n,ncol=p)
  curr.index = match(curr.data$curr,par.values)
  prev.index = match(curr.data$prev,par.values)
  for(i in 1:n){ X[i,curr.index[i]] = 1; X[i,prev.index[i]] = -1 }; 
  # Pseudo-outcome vector (n x 1)
  kappa = (y-0.5)
  
  #AR(1) Precision Matrix (p x p), alpha \propto Normal(0,tau Q)
  #We could use this as the initialization for Q#
  if(Qtype==1){
    Q = matrix(0,p,p); Q[1,1] = Q[p,p] = 1
    for(t in 2:(p-1)){ Q[t,t] = 1+0.99^2; Q[t,t+1] = Q[t,t-1] = -0.99; }
    Q[p,p-1] <- Q[1,2] <- -0.99 #we initialize phi = 0.99 in this case above#
    phi <- 0.99}else{ 
      Q = iCARQ(p); #specifies for the use of improper CAR of univariate dimension alpha
      Q = Q + diag(p)*lambda} #as suggested by Rue and Held page 108 "Weak Norm"
  
  # Scaling down the prior distribution for alpha by scaling the prior density for tau, by 
  # restricting the probability of the range of the alphas being less than 1 = 5%
  # more information provided in the paper. 
  
  samps = inla.qsample(n=10000,Q=Q,mu=rep(0,p),constr=list(A = matrix(rep(1, p), 1, p), e = 0)) #sample alpha with constraint
  numrang <- apply(samps, 2, function(x) max(x) - min(x)); rm(samps)
  x <- quantile(numrang, probs = .05)
  
  # Storage object
  alpha.samps = matrix(NA,nrow=nsamps,ncol=p)
  tau.samps = rep(NA,nsamps)
  if(Qtype == 2){
    phi.samps = rep(NA,nsamps) # Phi's being collected #
  }else{phi.samps = NA}
  # Initial values
  tau = 0.5 
  alpha = cbind(rnorm(p,0,1/sqrt(tau)))
  alpha.star = alpha - mean(alpha)
  
  
  # Gibbs Sampler  
  for(samp in 1:(nsamps+warmup)){
    # Sample alpha
    P0 = tau*Q
    omega = rpg.devroye(n,1,X%*%alpha)        # Update auxillary parameters

    V = chol2inv(chol(t(X)%*%(X*omega) + P0)) # Variance matrix, V has row and column sums = 0
    m = V%*%(t(X)%*%kappa)                    # Mean vector
    alpha = m + t(chol(V))%*%rnorm(ncol(X))#}  # Unconstrained
    
    #sample from posterior for tau
    tau = rgamma(1, shape = (a + p/2), rate = ((t(alpha)%*%Q%*%alpha)/2 + b/x^2)) 
    
    #Sample phi
    if(Qtype == 2){
      phi = Qreject(phi = phi, alpha = alpha, tau = tau)
      Q <- pluginQ(p, phi)
    }
    # Save posterior sample 
    if(samp>warmup){
      alpha.samps[samp-warmup,] = alpha-rowSums(V)*sum(alpha)/sum(V) #Corrected for Sum to Zero Constraint
      tau.samps[samp-warmup] = tau
      if(Qtype == 2){
        phi.samps[samp-warmup] = phi}
    } 
    
  }
  # Brier score evaluates how close the estimated probabilities of preferring the configurations
  # are to the true preferences given in the schedule (the observed Y's), for curr.data
  Brier <- BrierScore1D(alpha.samps, curr.data, par.values)
  # Return posterior samples
  return(list(alpha.samps=alpha.samps,tau.samps=tau.samps, phi.samps = phi.samps, Brier = Brier)) 
}



expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  
  do.call(rbind, lapply(seq_along(x), g))
}

# computes inverse-logits in differences between two preferences X and Y
expits <- function(x,y){val = exp(x - y)/(1+exp(x - y)); return(val)}

# Computes the RMSE between true and estimated preferences over the configurations
#fit: is the output from the model fit (from bicar2D or btgmrf)
#true.alpha: a vector of length equal to number of configurations, of true latent preferences
get.errors <- function(fit, true.alpha){
  nalpha <- ncol(fit$alpha.samps); alpha.post <- fit$alpha.samps
  crossed <- data.frame(expand.grid.unique(1:nalpha, 1:nalpha, include.equals = FALSE))
  total <- c(); totaltrue <- c(); 
  
  for(i in 1:(nalpha-1)){
    guy <- alpha.post[,i]; comps <- crossed[crossed$X1 == i,2]
    oneexp <- colMeans(sapply(comps, function(x)expits(guy, alpha.post[,x])))
    total <- c(total, oneexp);
    guy2 <- true.alpha[i]; 
    twoexp <- sapply(comps, function(x)expits(guy2, true.alpha[x]))
    totaltrue <- c(totaltrue, twoexp)
  }
  result <- sqrt(sum((total - totaltrue)^2)/choose(nalpha,2))
  return(list("RMSE" = result))
}

# Function that provides the configuration with the maximum posterior mean in preference, 
# and the expit between the current estimated best and the preference of the configuration with 
# the truly highest preference; 
#fit: is the output from the model fit (from bicar2D or btgmrf)
#true.alpha: a vector of length equal to number of configurations, of true latent preferences
get.BEST <- function(fit, true.alpha){
  nalpha <- ncol(fit$alpha.samps)
  bestX <- (1:nalpha)[which.max(colMeans(fit$alpha.samps))]
  bestdiff <- expits(true.alpha[bestX],max(true.alpha))
  return(list("AlphHatBest" = bestX, "BestPrefDiff" = bestdiff, "TrueBest" = which.max(true.alpha)))
}
#.5 will be the best, intuitively higher is better#



# Function that provides the Trial Quality metric, the average probability of 
# preferring the configurations tested over the true optimal configuration, over
# the entire device calibration history. 

#true.alpha:a vector of length equal to number of configurations, of true latent preferences
#curr.data: the accumulated monthly data. One column are the binary outcomes Y's for the patient-
# reported preferences between the consecutively tested configurations, provided in the next two 
# columns "curr" and "prev". 
# an example of curr.data will be provided 
#estim.alpha: alpha object from model fit (from bicar2D or btgmrf)
# usually this can be called if you named your model fit "fit" and 
# calling fit$alpha
# time: number of months intended (before trial stopping) - usually set to 12
# days: number of assumed days per month (we assumed 30)
get.quality <- function(true.alpha, curr.data, estim.alpha, time, days){
  hypothetical.time <- time * (days) #months times number of comparisons
  #need to include the first setting of each month, too. 
  observ.settings <- numeric(hypothetical.time); obs.days.tot <- length(curr.data[,3]); nmonths <- obs.days.tot/29 
  # see if this is below working
  for(j in 0:(nmonths-1)){observ.settings[29*j + j+1] <- curr.data[(29*j + 1),4]; observ.settings[29*j + j+1 + 1:29] <- curr.data[(29*j) + 1:29,3]}
  nalpha <- length(true.alpha); topsetting <- which.max(colMeans(estim.alpha))
  observ.settings[1:obs.days.tot] <- curr.data[,3]; observ.settings[(obs.days.tot+1):hypothetical.time] <- topsetting
  return(sum(expits(true.alpha[(observ.settings)],max(true.alpha)))/length(observ.settings))
}

# One simulation of a year long trial #
#par.values: vector of configuration indices (i.e., c(1:100))
#true.alpha: vector of true preferences for the configurations in the 1-D context
# b.size: how large would you like the test batch of configurations each month? (we set it to 8)
# b.retain: how many configurations would you like to keep from the previous month, 
# to be tested in the following month? (we set it to 2)
# days: number of days per month (we set this to 30)
# Qtype: if equal to 1, this sets phi = .99, if equal to 2, use random walk metropolis hastings
# to sample for phi. 
#lambda: a hack for invertible Q
#type: type of Batch Acquisition strategy to employ
# choose from: "EI", "Slice", "qEI", or "None"
#sampling: type of sampling mechanism from the acquisition function 
# if you chose "EI" or "Slice", you may choose "Penaldist" or "Slice"
# if you chose "qEI" specify "None" for sampling
# if you chose "None" for type, then specify "Thompson" for sampling 
# information about these types of Batch Acquisition strategies can be found
# in the article associated with this code. 

onesim <- function(par.values, true.alpha, b.size, b.retain, days, Qtype=1, lambda = 0.0001, type, sampling){
  time = 12; nsamps = 2000; burns = 500
  curr.data = NULL; BrierTraject <- matrix(0, ncol = 2, nrow = time);
  month = 1; collects = matrix(NA, ncol = 5, nrow = time)
  futmat <- matrix(NA, ncol = time, nrow = 1); supmat <- matrix(NA, ncol = time, nrow = 1)
  for(month in 1:time){
    #Get next batch and allocation order for next month
    if(month==1){batch = as.vector(quantile(par.values,0:7/7,type=1))} 
    if(month >1){batchfunc = acquisChoice1D(type, sampling, alpha.samps = fit$alpha.samps,
                                            par.values, month.past = month.data, b.retain, b.size);
    batch <- batchfunc[[1]]}
    order = get.order(nvalues=b.size,ndays=days)
    month.data = get.month.data2(alphatrue = true.alpha, batch=batch,order=order) 
    if(month > 1){Brier1 <- BrierScore1D(fit$alpha.samps, month.data,par.values)}
    #print(sprintf("Brier of Type 1 is equal to %f", Brier1))}
    if(month == 1){Brier1 <- 0}
    curr.data = rbind(curr.data,month.data)
    obs.points <- unique(curr.data$curr)
    fit = btgmrfADAM(par.values, curr.data=curr.data, Qtype = Qtype, lambda = lambda, nsamps = nsamps, warmup = burns)
    BrierTraject[month,] <- c(Brier1, fit$Brier)
    
    errors = get.errors(fit, true.alpha); RMSE = errors[[1]]; 
    best = get.BEST(fit, true.alpha); estbest = best[[1]]; offness= best[[2]]
    burden <- get.quality(true.alpha, curr.data, fit$alpha.samps, month, days)
    collects[month,] <- c(estbest, month, RMSE, offness, burden) 
    futmat[month] <- Futility(fit$alpha.samps)
    supmat[month] <- Superiority(fit$alpha.samps)
  }
  return(list("SimMeasures" = collects,"BrierTraject" = BrierTraject,
              "FutilityMat" = futmat, "SuperiorMat" = supmat))
}

# One simulation of a year long trial in the 2-D context #
#x1.values: vector of indices for levels of the first parameter of the device (i.e. c(1:10))
#x2.values: vector of indices for levels of the second parameter of the device (i.e. c(1:11))
#true.alpha: vector of true preferences for the configurations in the 2-D context (in the examples
# previously provided for x1.values and x2.values, the true.alpha vector would be of length 110, 
#i.e, 10 x 11)

# b.size: how large would you like the test batch of configurations each month? (we set it to 8)
# b.retain: how many configurations would you like to keep from the previous month, 
# to be tested in the following month? (we set it to 2)
# days: number of days per month (we set this to 30)
#model.type: 1 specifies the 2NCAR model while 2 specifies the additive 2-D model
#type: type of Batch Acquisition strategy to employ
# choose from: "EI", "Slice", "qEI", or "None"
#sampling: type of sampling mechanism from the acquisition function 
# if you chose "EI" or "Slice", you may choose "Penaldist" or "Slice"
# if you chose "qEI" specify "None" for sampling
# if you chose "None" for type, then specify "Thompson" for sampling 
# information about these types of Batch Acquisition strategies can be found
# in the article associated with this code. 

onesim2D <- function(x1.values, x2.values, true.alpha, b.size, b.retain, days, model.type = 1, type, sampling){
  nsamps = 2000; burns = 500
  curr.data = c(); flag = 0; month = 1; time = 12; collects <- matrix(NA, nrow = time, ncol = 5)
  BrierTraject <- matrix(0, ncol = 2, nrow = time);
  J = length(x1.values); K = length(x2.values); par.values = 1:(J*K)
  futmat <- matrix(NA, ncol = time, nrow = 1); supmat <- matrix(NA, ncol = time, nrow = 1)
  for(month in 1:time){
    #Get next batch and allocation order for next month
    if(month==1){batch = as.vector(quantile(par.values,0:7/7,type=1))} 
    if(month >1){batchfunc = acquisChoice2D(type, sampling, alpha.samps = fit$alpha.samps,
                                            x1.values, x2.values, month.past = month.data, b.retain = 2, b.size =8);
    batch <- batchfunc[[1]]}
    order = get.order(nvalues=b.size, ndays=days)
    month.data = get.month.data2D(alphatrue = true.alpha, batch=batch,order=order); 
    month.data <- Data1Dto2D(month.data, x1.values, x2.values)
    if(month > 1){Brier1 <- BrierScore(fit$alpha.samps, month.data, x1.values, x2.values)}
    if(month == 1){Brier1 <- 0}
    curr.data = rbind(curr.data,month.data)
    obs.points <- unique(curr.data[,c(2,3)])
    #print(sprintf("Commencing bi-CAR model for Month %s", month))
    if(model.type == 1){
      fit = biCAR2D(curr.data, x1.values, x2.values, nsamps = 2000, warmup = 500)}
    if(model.type == 2){
      #sprintf(print("Initiating Additive 2-D Model"))
      fit = fit.additive2D.model(dat = curr.data, x1.values, x2.values,nsamps=2000,warmup=burns)
    }
    BrierTraject[month,] <- c(Brier1, fit$Brier)
    
    errors = get.errors(fit, true.alpha); RMSE = errors[[1]]
    best = get.BEST(fit, true.alpha); estbest = best[[1]]; offness= best[[2]]
    burden <- get.quality(true.alpha, curr.data, fit$alpha.samps, month, days)
    collects[month,] <- c(estbest, month, RMSE, offness, burden) 
    futmat[month] <- Futility(fit$alpha.samps); 
    supmat[month] <- Superiority(fit$alpha.samps)
  }
  #print(sprintf("Simulation %s Done", i))
  #}
  return(list("SimMeasures" = collects,"BrierTraject" = BrierTraject,
              "FutilityMat" = futmat, "SuperiorMat" = supmat))
}

#Function to fit the 2NCAR model for binary outcome#
#dat: curr.dat; a matrix with the a $y column, the binary indicators of patient reported preference
# and the next two columns labeled $x1.prev and $x2.prev, where these correspond to the 
# 2-D configuration tested on the previous day, and the last two columns corresponding to the 
# the 2-D configuration tested on the current day ($x1.curr, $x2.curr)
#x1.values: vector of indices for levels of the first parameter of the device (i.e. c(1:10))
#x2.values: vector of indices for levels of the second parameter of the device (i.e. c(1:11))
library(BayesLogit)
biCAR2D = function(dat, x1.values, x2.values, nsamps=2000,warmup=500,
                   hyp = c(5/2, 5/2)){
  
  J = length(x1.values); K = length(x2.values)
  n = nrow(dat); y = dat$y; x1.curr = dat$x1.curr; x2.curr = dat$x2.curr;
  x1.prev = dat$x1.prev; x2.prev = dat$x2.prev
  
  p = J*K
  a = hyp[1]; b = hyp[2]; 
  # Design matrix (n x JK)
  X = matrix(0,nrow=n,ncol=J*K)
  for(i in 1:n){
    X[i,x2.curr[i]+(x1.curr[i]-1)*K] = 1
    X[i,x2.prev[i]+(x1.prev[i]-1)*K] = -1
  } 
  # Pseudo-outcome vector (n x 1)
  kappa = cbind(y-0.5)
  
  #AR(1) Precision Matrix (p x p), alpha \propto Normal(0,tau Q)
  
  QJ = matrix(0,J,J); diag(QJ) = c(1,rep(2,J-2),1)
  QJ[cbind(2:J,2:J-1)] = QJ[cbind(2:J-1,2:J)] = -1
  QK = matrix(0,K,K); diag(QK) = c(1,rep(2,K-2),1)
  QK[cbind(2:K,2:K-1)] = QK[cbind(2:K-1,2:K)] = -1
  
  # Storage object
  alpha.samps = matrix(NA,nrow=nsamps,ncol=p)
  tau.samps = rep(NA,nsamps)
  phi.samps = rep(NA,nsamps) # Phi's being collected #
  # Initial values
  tau = 0.5 
  alpha.vec = rnorm(p,0,1)
  A = QK %x% diag(J); B = diag(K) %x% QJ; phi = .5
  Q0 = 2*(1-phi)*(A)+2*(phi)*(B) 
  #b0 <- numeric(J*K); b0[1] <- 1; b0[J*K] <- -1
  diag(Q0) = diag(Q0) +  1e-10 #Prior Precision sans tau
  accept <- 0
  
  samps = inla.qsample(n=10000,Q=Q0,mu=rep(0,J*K),constr=list(A = matrix(rep(1, J*K), 1, J*K), e = 0)) #sample alpha with constraint
  numrang <- apply(samps, 2, function(x) max(x) - min(x)); rm(samps)
  x <- quantile(numrang, probs = .05)
  
  for(samp in 1:(nsamps+warmup)){
    # Sample alpha
    omega = rpg.devroye(n,1,X%*%alpha.vec) # update auxillary parameters
    Q0 = 2*phi*(B)+2*(1-phi)*(A) 
    diag(Q0) = diag(Q0) + + 1e-10 #Prior Precision sans tau
    Q = t(X)%*%(X*omega)+tau*Q0;
    Q = inla.as.sparse(Q)
    #m = inla.qsolve(Q,diag(J*K),method="solve")%*%(t(X)%*%kappa) # Posterior Precision and Mean
    m = inla.qsolve(Q,diag(J*K),method="solve")%*%(crossprod(X,kappa)) # Posterior Precision and Mean
    alpha.vec = inla.qsample(n=1,Q=Q,mu=m,constr=list(A = matrix(rep(1, J*K), 1, J*K), e = 0)) #sample alpha with constraint
    
    #wrapped Q0 multiplication to speed things up, maybe... 
    tau = rgamma(1, shape = (a + p/2-1), rate = c((t(alpha.vec)%*%(Q0%*%alpha.vec))/2 + b/x^2)) #if we use iCAR modeling scheme
    
    rej.stuff <- Reject.Psi(phi, A, B, J, K, alpha.vec, tau); phi <- rej.stuff[[1]]; accept <- accept + rej.stuff[[2]]
    
    
    if(samp>warmup){
      #if(samp%%100 ==0 ){print(sprintf("Sample Collected %d", samp-warmup))}
      #if((samp-warmup) %in% quantile(1:nsamps,1:9/10,type=1)) print(paste("Posterior Sampling ",round(100*(samp-warmup)/nsamps),"% Complete",sep=""))
      alpha.samps[samp-warmup,] = alpha.vec
      tau.samps[samp-warmup] = tau
      phi.samps[samp-warmup] = phi
    }
  } 
  Brier <- BrierScore(alpha.samps, dat, x1.values, x2.values)
  # Return posterior samples
  return(list(alpha.samps=alpha.samps,tau.samps=tau.samps, phi.samps = phi.samps, Brier = Brier, accept = accept/nsamps)) 
}

#Have a rejection sampler with sort of a random walk (may change later)
#The random walk adds a random uniform(0,1) to the previous value for phi#
Qreject2D = function(phi, v, alpha, tau, Eig1, Eig2, p1, p2){
  propPhi = runif(1) #random walk
  #safe guarding against problematic values of phi
  if(logLike4phi(tau, phi = propPhi, phi2 = phi, v, alpha, Eig1, Eig2, p1, p2) - logLike4phi(tau, phi = phi, phi2 = propPhi, v, alpha, Eig1, Eig2, p1, p2) > log(runif(1))){
    newPhi <- propPhi; accept = 1}else{newPhi <- phi; accept = 0}
  return(list("Phi" = newPhi, "Accept" = accept))
}

# This function fits the additive 2D model, the alternative model to the 2NCAR 
#dat: curr.dat; a matrix with the a $y column, the binary indicators of patient reported preference
# and the next two columns labeled $x1.prev and $x2.prev, where these correspond to the 
# 2-D configuration tested on the previous day, and the last two columns corresponding to the 
# the 2-D configuration tested on the current day ($x1.curr, $x2.curr)
#x1.values: vector of indices for levels of the first parameter of the device (i.e. c(1:10))
#x2.values: vector of indices for levels of the second parameter of the device (i.e. c(1:11))
fit.additive2D.model = function(dat,x1.values,x2.values,nsamps=2000,warmup=500){
  
  #Data pre-processing
  n = nrow(dat); J = length(x1.values); K = length(x2.values)
  y = dat$y; x1.curr = dat$x1.curr; x2.curr = dat$x2.curr; x1.prev = dat$x1.prev; x2.prev = dat$x2.prev
  
  # Pseudo-outcome vector (n x 1)
  kappa = cbind(y-0.5)
  
  # Design matrix (n x [J+K])
  XJ = t(apply(cbind(x1.curr,x1.prev),1,function(x) 1*(x[1]==x1.values)-1*(x[2]==x1.values)))
  XK = t(apply(cbind(x2.curr,x2.prev),1,function(x) 1*(x[1]==x2.values)-1*(x[2]==x2.values)))
  X = cbind(XJ,XK)
  
  # Prior Precision Matrices (J x J) and (K x K)
  QJ = matrix(0,J,J); diag(QJ) = c(1,rep(2,J-2),1) + 1e-10; 
  QJ[cbind(2:J,2:J-1)] = QJ[cbind(2:J-1,2:J)] = -1
  QK = matrix(0,K,K); diag(QK) = c(1,rep(2,K-2),1) + 1e-10; 
  QK[cbind(2:K,2:K-1)] = QK[cbind(2:K-1,2:K)] = -1
  
  # Find prior distribution for tau for individual additive components 
  samps = inla.qsample(n=10000,Q=QK,mu=rep(0,K),constr=list(A = matrix(rep(1, K), 1, K), e = 0)) #sample alpha with constraint
  numrang <- apply(samps, 2, function(x) max(x) - min(x)); rm(samps)
  xK <- quantile(numrang, probs = .05)
  
  samps = inla.qsample(n=10000,Q=QJ,mu=rep(0,J),constr=list(A = matrix(rep(1, J), 1, J), e = 0)) #sample alpha with constraint
  numrang <- apply(samps, 2, function(x) max(x) - min(x)); rm(samps)
  xJ <- quantile(numrang, probs = .05)
  
  # Storage objects
  gamma.samps = matrix(NA,nrow=nsamps,ncol=J)
  lambda.samps = matrix(NA,nrow=nsamps,ncol=K)
  tau.samps = matrix(NA,nrow=nsamps,ncol=2)
  
  # Initial values
  tau = c(1,1)
  gamma.temp = cbind(rnorm(J)); gamma = gamma.temp - mean(gamma.temp)
  lambda.temp = cbind(rnorm(K)); lambda = lambda.temp - mean(lambda.temp)
  
  # Gibbs Sampler  
  for(samp in 1:(nsamps+warmup)){
    # Sample alpha
    P0 = rbind(cbind(tau[1]*QJ,matrix(0,J,K)),cbind(matrix(0,K,J),tau[2]*QK))
    omega = rpg.devroye(n,1,X%*%cbind(c(gamma,lambda))) # update auxillary parameters
    V = chol2inv(chol(t(X)%*%(X*omega) + P0)) # Variance matrix
    m = V%*%(t(X)%*%kappa)                    # Mean vector
    theta.temp = m + t(chol(V))%*%rnorm(ncol(X))
    gamma.temp = theta.temp[1:J,]; lambda.temp = theta.temp[J+1:K,] #Unconstrained
    gamma = gamma.temp-rowSums(V[1:J,1:J])*sum(gamma.temp)/sum(V[1:J,1:J]) #Corrected for Sum to Zero Constraint
    lambda = lambda.temp-rowSums(V[J+1:K,J+1:K])*sum(lambda.temp)/sum(V[J+1:K,J+1:K]) #Corrected for Sum to Zero Constraint
    
    # Sample taus
    tau.gamma = rgamma(1,shape=(5+J-1)/2,rate=((5/2)/xJ^2 + t(gamma)%*%QJ%*%gamma)/2)
    tau.lambda = rgamma(1,shape=(5+K-1)/2,rate=((5/2)/xK^2 +t(lambda)%*%QK%*%lambda)/2)
    
    # Save posterior sample 
    if(samp>warmup){
      gamma.samps[samp-warmup,] = gamma
      lambda.samps[samp-warmup,] = lambda
      tau.samps[samp-warmup,] = c(tau.gamma,tau.lambda)
    } 
  }
  # problem here --> ordering #
  alpha.samps <- (t(sapply(1:nsamps, function(m) as.vector(outer(gamma.samps[m,], lambda.samps[m,], FUN = "+")))))
  Brier <- BrierScore(alpha.samps, dat, x1.values, x2.values)
  # Return posterior samples
  return(list(gamma.samps=gamma.samps,lambda.samps=lambda.samps, alpha.samps = alpha.samps, tau.samps=tau.samps, Brier = Brier)) 
}

#changes the data from 1D (for 2ncar model) to 2D (for additive model)
Data1Dto2D <- function(dat1D, x1.values, x2.values){
  lattice <- referenceLattice(length(x1.values), length(x2.values))
  new2D <- matrix(0, ncol = 5, nrow = nrow(dat1D)); new2D[,1] <- dat1D$y;
  for(j in 1:nrow(new2D)){
    curr2 <- which(lattice==dat1D[j,3], arr.ind=TRUE)
    prev2 <- which(lattice==dat1D[j,4], arr.ind=TRUE)
    new2D[j,c(2:5)] <- cbind(curr2, prev2)
  }
  new2D <- as.data.frame(new2D)
  names(new2D) <- c("y", "x1.curr", "x2.curr", "x1.prev", "x2.prev")
  return(new2D)
}

#converts 1D to 2D configuration
oneDto2D <- function(oneDvals, x1.values, x2.values){
  lattice <- referenceLattice(length(x1.values), length(x2.values))
  new2D <- matrix(0, ncol = 2, nrow = length(oneDvals));
  for(j in 1:nrow(new2D)){
    new2D[j,] <- which(lattice==oneDvals[j], arr.ind=TRUE)
  }
  new2D <- as.data.frame(new2D)
  return(new2D)
}

#to convert 2 dim setting pairs into 1-D #
TwoDto1D <- function(TwoDat, x1.values, x2.values){
  refs <- referenceLattice(length(x1.values),length(x2.values))
  setts <- sapply(1:nrow(TwoDat), function(x) refs[TwoDat[x,1], TwoDat[x,2]])
  return(setts)
}

# Function for cutting the trial early for the
# Calibration Convergence Rule with margin of convergence of 0.5
Superiority <- function(alphasamps){
  alphastar <- which.max(colMeans(alphasamps))
  dimsd <- dim(alphasamps); long <- dimsd[1]
  Sup.Counts <- sapply(1:long, function(x){prob <- (max(alphasamps[x,]) - alphasamps[x,alphastar] < .5)})
  outt <- mean(Sup.Counts)
  return(outt)
}

#Function for cutting trial early for the
# Preference Neutrality Rule with range cut off of 1
Futility <- function(alpha.samps){
  dimsd <- dim(alpha.samps); long <- dimsd[1]
  Fut.Counts <- sapply(1:long, function(x){prob <- (max(alpha.samps[x,]) - min(alpha.samps[x,]) < 1)})
  outt <- mean(Fut.Counts);
  return(outt)
}

#Unnormalized posterior of phi parameter for 2NCAR
phi.unnorm.fc.log = function(phi,alpha,tau, J, K, A, B){
  eigsum = sum(unlist(sapply(apply(expand.grid(1:K,1:J),1,function(x) 4*(phi*tau*(1-cos(pi*(x[1]-1)/K))+(1-phi)*(1-cos(pi*(x[2]-1)/J)))),
                             function(x){y <- (log(x[x > 0]))})))  #the sum of all nonzero eigenvalues 
  P0 = 2*phi*(A)+2*(1-phi)*(B) + 1e-10*diag(J*K)
  d = (0.5*(eigsum-tau*as.numeric(t(alpha)%*%(P0%*%alpha))))
  return(as.numeric(d))
} 

#Unnormalized posterior of psi, the transformed phi
psi.unnorm.fc.log = function(psi, A, B, J, K, alpha, tau){
  phi = 1/(1 + exp(-psi))
  psi.log.dens <- phi.unnorm.fc.log(phi, alpha, tau, J, K, A,B) + log((exp(-psi)/(1 + exp(-psi))^(2)))
  return(psi.log.dens)
}

# Random Jitter for psi, transforms back to phi, metropolis hasting's step
Reject.Psi <- function(phi, A, B, J, K, alpha, tau){
  psi <- log(phi/ (1-phi))
  psi.prop = psi + rnorm(1, 0, 2); 
  if(psi.unnorm.fc.log(psi.prop, A, B, J, K, alpha, tau) - psi.unnorm.fc.log(psi, A, B, J, K, alpha, tau) > log(runif(1))){
    psi = psi.prop; accept = 1}else{accept = 0}
  phiout <- 1/(1+exp(-psi))
  return(list("Phi" = phiout, "Accept" = accept))
}

#vectorizes a 2D lattice indexes#
corrds <- function(J, j, k){ 
  vals <- c()
  for (i in 1:length(j)){
    vals <- c(vals, j + J*(k-1))
  }
  return(vals)
}

#BS1 is model fit to data from months 1,...,m-1 
#predicting data from month m. So you predict month 2 based on month 1 fit
#, month 3 based on month 1,2 fit, etc.

#BS2 is model fit to data from months 1,...,m predicting data from months
#1,...,m. So you predict month 1 based on month 1 fit, months 1&2 based on
#months 1&2 fit, etc.

#Brier Score evaluations for 1-D
BrierScore1D <- function(alphasamps, monthdata, par.values){
  # If you input months 1,...m-1, then you are predicting responses of month m,
  # use fit(curr.data up till m-1) to predict month.data month m's (y)
  # this is dependent on different inputs
  # alphasamps are the alpha posterior samples based on months 1,...,m-1, and monthdata should be for
  # month m
  
  #If you want the second type of Brier score, input posterior samples for alpha evaluated 
  # from accrued data for months 1,...m, and then place ALL month data (the stacked monthly data)
  # up till month m. 
  
  BrierType <- mean((monthdata$y - mapply(function(x,y) expits(x,y), colMeans(alphasamps)[monthdata$curr],colMeans(alphasamps)[monthdata$prev]))^2)
  return(BrierType)
}

#Brier Score evaluations 2-D
BrierScore <- function(alphasamps, monthdat, x1.values, x2.values){
  # If you input months 1,...m-1, then you are predicting responses of month m,
  # use fit(curr.data up till m-1) to predict month.data month m's (y)
  # this is dependent on different inputs
  # alphasamps are the alpha posterior samples based on months 1,...,m-1, and monthdata should be for
  # month m
  
  #If you want the second type of Brier score, input posterior samples for alpha evaluated 
  # from accrued data for months 1,...m, and then place ALL month data (the stacked monthly data)
  # up till month m. 
  
  reflat = referenceLattice(length(x1.values), length(x2.values))
  Dim1 <- Convert2Dto1D(monthdat, reflat)
  BrierType <- mean((Dim1$y - mapply(function(x,y) expits(x,y), colMeans(alphasamps)[Dim1$curr],colMeans(alphasamps)[Dim1$prev]))^2)
  return(BrierType)
}

# Function to sample 1-Mode 2-D latent preference surfaces 
#J: number of values for parameter 1
#K: number of values for parameter 2 of the device
#mu: vector of length 2 which resides in the [-3,-3] X [-3,3] X [3,3] X [3,-3] grid
# that reflects the mean of the 2-D gaussian distribution location 
# Sigma: 2x2 covariance matrix, make sure it is positive definite 
sample2DGaussUni <- function(J, K, mu, Sigma){
  stor <- matrix(0, ncol = J, nrow = K)
  for(j in 1:J){
    for (k in 1:K){
      stor[k, j] <- (dmvnorm(6*c(j/J - .5, k/K - .5), mu, Sigma))  
    }
  }
  stor <- 40*stor - 2
  return(stor)
}

########################### Batch Bayesian Optimization Functions ##########################
#######################Distance function for penalization ##################################

#Euclidian Distance function for 2-D
#dat 1 needs to be a N x 3 matrix of 2-D configurations, with last column being estimated mean preference
# for the corresponding configuration
#dat 2 needs to be a N x 3 matrix of 2-D configurations, with last column being estimated mean preference
# for the corresponding configuration 
#l2: if dat2 is a column vector, then specify l2 = 1
# dat 1 is taken to be the "candidate values, i.e., not tested yet" and dat2 is the batch elements for month m
distfunc2D <- function(dat1, dat2, l2 = 0){
  if(l2 == 0){
    dists <- sapply(1:nrow(dat1), function(x){sapply(1:nrow(dat2), function(l){sqrt(sum((dat1[x,1] - dat2[l,1])^2 + (dat1[x,2]-dat2[l,2])^2 + (dat1[x,3]-dat2[l,3])^2 ))})})}
  if(l2 == 1){
    dists <- sapply(1:nrow(dat1), function(x){sqrt(sum((dat1[x,1] - dat2[1])^2 + (dat1[x,2]-dat2[2])^2 + (dat1[x,3]-dat2[3])^2 ))})}
  return(dists)
}


#Euclidian Distance function for 1-D
#dat 1 needs to be a N x 2 matrix of 1-D configurations, with last column being estimated mean preference
# for the corresponding configuration
#dat 2 needs to be a N x 2 matrix of 1-D configurations, with last column being estimated mean preference
# for the corresponding configuration 
#l2: if dat2 is a column vector, then specify l2 = 1
# dat 1 is taken to be the "candidate values, i.e., not tested yet" and dat2 is the batch elements for month m

distfunc1D <- function(dat1, dat2, l2 = 0){
  if(l2 == 0){
    dists <- sapply(1:nrow(dat1), function(x){sapply(1:nrow(dat2), function(l){sqrt(sum((dat1[x,1] - dat2[l,1])^2 + (dat1[x,2]-dat2[l,2])^2))})})
  }
  if(l2 == 1){
    dists <- sapply(1:nrow(dat1), function(x){sqrt(sum((dat1[x,1] - dat2[1])^2 + (dat1[x,2]-dat2[2])^2))})
  }
  return(dists)
}


###################################### Batch Acquisition Functions ##########################
# 2D - Generalized Slice Sampler with Inverse Distance Penalization

#acquis: matrix of size J x K evaluated acquisitions (through EI or UCB)
# x1.values: 1:J
# x2.values: 1:K
# month.past: past month's schedule and binary outcomes  
# b.retain: how many configurations do you want to keep from the previous month's batch
# b.size: batch size  (we set this to 8)
# toppast: is the previous month m-1's configuration INDEX 
# with the maximum posterior mean in latent preference
gSlice2DPenalizDist <- function(acquis, x1.values, x2.values, month.past, b.retain = 2, b.size = 8, toppast){
  
  b.grab <- b.size - b.retain
  acquis[acquis < 0] <- log(1 + exp(acquis[acquis < 0])) #force acquisition values to be non-negative
  grid.values <- expand.grid(x1.values, x2.values) 
  J <- max(grid.values[,1]); K <- max(grid.values[,2])
  longacq <- c(acquis) #vectorize the matrix of acquisitions
  info <- cbind(expand.grid(1:J, 1:K), longacq); uniquepast <- unique(month.past[,c(2,3)])
  pastinfo <- which(rowSums(apply(uniquepast[,c(1,2)], 1, function(a) apply(info[,c(1,2)], 1, function(b) all(a==b)))) == 1) #converts it into 1-D setting values
  
  # Obtaining the past month's past top configuration 
  # to place it into the batch 
  if(!is.null(toppast)){
    pastacqui <- max(longacq[pastinfo][pastinfo != toppast])
    S <- matrix(NA, nrow = b.size, ncol= 3); S[1,] <- as.numeric(info[toppast,])
  }
  if(is.null(toppast)){
    toppast <- which.max(longacq[pastinfo]);  pastacqui <- max(longacq[pastinfo][pastinfo != toppast])
    S <- matrix(NA, nrow = b.size, ncol= 3); S[1,] <- as.numeric(info[toppast,])
  }
  
  #Filling in the batch 
  for(i in 2:b.size){
    # from previous month's batch
    if(i <= b.retain & i > 1){ 
      
      dist <- distfunc2D(info[pastinfo,], S[i-1,], l2 = 1) #evaluate 2D distance function between 
      # the past month's batch elements and the previous batch elements 
      
      dist[dist == 0] <- 1; dist <- 1 - (1/dist) #if the distance is 0, then no penalty
      preacqui <- longacq[pastinfo]*dist # apply the penalty 
      distmat <- matrix(1, nrow = J, ncol = K); 
      for(j in 1:length(pastinfo)){coordsd <- c(info[pastinfo[j],1], info[pastinfo[j],2]); distmat[coordsd[1], coordsd[2]] <- dist[j]}
      acquis <- acquis * distmat # element wise multiplication of the penalties on all acquisitions, 
      # not just previous month's
      u1 <- runif(1, min(preacqui), max(preacqui)) #prefix the bottom to be 0; can't sample the last 
      #acquisition value evaluated at 0-valued previous best setting option so
      # just sample from maximum of the previous months' discredited values 
      accept = 0; while(accept == 0){ #mini sample from past month data
        # slice sampling (see Neal 2003)
        s0 <- as.numeric(info[sample(pastinfo, size = 1, replace= F),c(1,2)])
        alphaprop <- acquis[s0[1], s0[2]]
        if(alphaprop > u1){S[i,] <- c(as.numeric(s0[1]), as.numeric(s0[2]), alphaprop);
        accept = 1; pastacqui <- alphaprop}else{accept = 0; u1 <- u1*.85}
      }
    }
    zeroC <- 1
    # Open up to entire configuration space once you pass b.retain batch elements 
    if(i > b.retain){
      longacq <- c(acquis) #first penalty is in place on acquis
      info <- cbind(grid.values, longacq) 
      dist <- distfunc2D(info, S[i-1,], l2 = 1); #distance formula for (x,y,alpha) between past batch element
      #and the rest of the data including their past penalties 
      dist[dist == 0] <- 1; dist <- 1 - 1/dist #farther away implies less penalty
      acquis <- acquis * matrix(dist, ncol = K, nrow = J, byrow = F) #sequential punish
      u1 <- runif(1, 0, max(acquis))
      accept = 0; while(accept == 0){
        si <- grid.values[sample(1:nrow(grid.values), size = 1, replace= F),] #sample coords
        alphaprop <- acquis[as.numeric(si[1]),as.numeric(si[2])]# see last line of this code
        if(alphaprop > u1){S[i,] <- c(as.numeric(si[1]), as.numeric(si[2]), alphaprop);
        accept = 1}else{accept = 0; u1 <- u1*.85;zeroC <- zeroC + 1}
        if(zeroC > 150){S[i,] <- c(as.numeric(grid.values[sample(1:nrow(grid.values), size = 1, replace= F),]),.00001); break}
      }
    }
  }
  return(list("NewBatch" = S[,c(1,2)])) #return the batch from generalized penalized distance acquisition
}


# 1D - Generalized Slice Sampler with Inverse Distance Penalization

#acquis: matrix of size J x K evaluated acquisitions (through EI or UCB)
# par.values: configuration index, (i.e., 1:100)
# month.past: past month's schedule and binary outcomes  
# b.retain: how many configurations do you want to keep from the previous month's batch
# b.size: batch size  (we set this to 8)
# toppast: is the previous month m-1's configuration INDEX 
# with the maximum posterior mean in latent preference

#see annotations for gSlice2DPenalizDist for more detailed annotations
gSlice1DPenalizDist <- function(acquis, par.values, month.past, b.retain = 2, b.size = 8, toppast){
  
  b.grab <- b.size - b.retain; 
  acquis[acquis < 0] <- log(1 + exp(acquis[acquis < 0]))
  J = length(par.values)
  info <- cbind(par.values, acquis); uniquepast <- unique(month.past[,c(3)])
  #slice sample
  if(!is.null(toppast)){ #in the case of using EI surface 
    pastacqui <- max(acquis[uniquepast][uniquepast != toppast])
    S <- matrix(NA, nrow = b.size, ncol= 2); S[1,] <- c(toppast, acquis[toppast])}
  if(is.null(toppast)){ #in the case of using UCB acquisition surface
    toppast <- which.max(acquis[uniquepast]); pastacqui <- max(acquis[uniquepast][uniquepast != toppast])
    S <- matrix(NA, nrow = b.size, ncol= 2); S[1,] <- c(toppast, acquis[toppast])
  }
  exclude <- c()
  minacqui = 0
  for(i in 2:b.size){
    if(i <= b.retain & i > 1){ # b.retain # element should be determined from only previous month info
      
      dist <- distfunc1D(info[uniquepast,], S[i-1,], l2 = 1)
      dist[dist == 0] <- 1; dist <- 1 - (1/dist)
      preacqui <- acquis[uniquepast]*dist #previous batch discounted values from first chosen setting
      distvec <- rep(1, J) ; pastinfo <- cbind(info[uniquepast,], dist)
      for(j in 1:nrow(pastinfo)){distvec[pastinfo[j,1]] <- dist[j]}; 
      acquis <- acquis * distvec
      u1 <- runif(1, min(preacqui), max(preacqui)) #prefix the bottom to be 0; can't sample the last 
      #acquisition value evaluated at 0-valued previous best setting option so
      # just sample from maximum of the previous months' discredited values 
      accept = 0; while(accept == 0){ #mini sample from past month data
        s0 <- sample(pastinfo[,1], size = 1, replace= F)
        alphaprop <- acquis[s0]
        if(alphaprop > u1){S[i,] <- c(s0, alphaprop);
        accept = 1; pastacqui <- alphaprop}else{accept = 0; u1 <- u1*.85}
      }
    }
    #first penalty is in place on acquis
    #penalize by inverse distance
    if(i > b.retain){
      info <- cbind(par.values, acquis)
      dist <- distfunc1D(info, S[i-1,], l2 = 1); #distance formula for (x,y,alpha) between past batch element
      #and the rest of the data including their past penalties 
      dist[dist == 0] <- 1; dist <- 1 - 1/dist #farther away implies less penalty
      dist[S[1:(i-1),1]] <- 0
      acquis <- acquis * dist #sequential punish
      u1 <- runif(1, 0, max(acquis))
      accept = 0; while(accept == 0){
        si <- sample(par.values, size = 1, replace= F) #sample coords
        alphaprop <- acquis[si]# see last line of this code
        if(alphaprop > u1){S[i,] <- c(si, alphaprop);
        accept = 1}else{accept = 0; u1 <- u1*.85}}
      #below we just removed the past values, but now we can punish them
      #grid.values <- grid.values[-which(rownames(grid.values) == rownames(si)),]
      #pastacqui <- acquis[as.numeric(si[1]), as.numeric(si[2])] <- sometimes cannot beat this
      #print(sprintf("batch element %s found", i))
    }
  }
  return(list("NewBatch" = S[,c(1,2)]))
}

# Generalized Slice Sampler without Inverse Distance Penalization 

#acquis: matrix of size J x K evaluated acquisitions (through EI or UCB)
# x1.values: 1:J
# x2.values: 1:K
# month.past: past month's schedule and binary outcomes  
# b.retain: how many configurations do you want to keep from the previous month's batch
# b.size: batch size  (we set this to 8)
# toppast: is the previous month m-1's configuration INDEX 
# with the maximum posterior mean in latent preference

# More detailed annotations for this function are included in gSlice2DPenalizDist
# ignoring the sequential penalization 
# This function, instead of penalizing, simply removes previously chosen batch items from
# the candidate set 
gSlice2D <- function(acquis, x1.values, x2.values, month.past, b.retain = 2, b.size = 8, toppast){
  
  b.grab <- b.size - b.retain
  acquis[acquis < 0] <- log(1 + exp(acquis[acquis < 0]))
  grid.values <- expand.grid(x1.values, x2.values)  
  J <- max(grid.values[,1]); K <- max(grid.values[,2])
  longacq <- c(acquis)
  info <- cbind(expand.grid(1:J, 1:K), longacq); uniquepast <- unique(month.past[,c(2,3)])
  pastinfo <- which(rowSums(apply(uniquepast[,c(1,2)], 1, function(a) apply(info[,c(1,2)], 1, function(b) all(a==b)))) == 1) #converts it into 1-D setting values
  if(!is.null(toppast)){
    pastacqui <- max(longacq[pastinfo][pastinfo != toppast])
    S <- matrix(NA, nrow = b.size, ncol= 3); S[1,] <- as.numeric(info[toppast,])
  }
  if(is.null(toppast)){
    toppast <- pastinfo[which.max(longacq[pastinfo])];  pastacqui <- max(longacq[pastinfo][pastinfo != toppast])
    S <- matrix(NA, nrow = b.size, ncol= 3); S[1,] <- as.numeric(info[toppast,])
  }
  exclude = c(); 
  for(i in 2:b.size){
    if(i <= b.retain & i > 1){
      tt <- as.numeric(proc.time()[1])
      # b.retain # element should be determined from only previous month info
      u1 <- runif(1, 0, pastacqui) #prefix the bottom to be 0; can't sample the last 
      #acquisition value evaluated at 0-valued previous best setting option so
      # just sample from maximum of the previous months' discredited values 
      accept = 0; while(accept == 0){ #mini sample from past month data
        toppast <- c(toppast, exclude)
        s0i <- (info[sample(pastinfo[!(pastinfo %in% exclude)], size = 1, replace= F),c(1,2)])
        s0 <- as.numeric(s0i)
        alphaprop <- acquis[s0[1], s0[2]]
        if(alphaprop > u1){S[i,] <- c(as.numeric(s0[1]), as.numeric(s0[2]), alphaprop);
        accept = 1; exclude <- as.numeric(rownames(s0i)); 
        grid.values <- grid.values[-which(rownames(grid.values) == exclude),]}else{accept = 0; u1 <- u1*.85}
        pastacqui <- alphaprop;
        if(tt - as.numeric(proc.time()[1]) > 30){
          break
        }}
      
    }
    zeroC <- 1
    if(i > b.retain){
      tt <- as.numeric(proc.time()[1])
      u1 <- runif(1, 0, pastacqui);
      accept = 0; while(accept == 0){
        si <- grid.values[sample(1:nrow(grid.values), size = 1, replace= F),]
        alphaprop <- acquis[as.numeric(si[1]),as.numeric(si[2])]
        if(alphaprop > u1){S[i,] <- c(as.numeric(si[1]), as.numeric(si[2]), alphaprop);
        accept = 1}else{accept = 0; u1 <- u1*.85; zeroC <- zeroC + 1}
        if(zeroC > 150){S[i,c(1,2)] <- as.numeric(grid.values[sample(1:nrow(grid.values), size = 1, replace= F),]); break
        } }
      grid.values <- grid.values[-which(rownames(grid.values) == rownames(si)),]
      #print(sprintf("Collection of Item %s",i))
      pastacqui <- acquis[as.numeric(si[1]), as.numeric(si[2])];}
    if(sum(!is.na(S[,1])) == b.size){break}
  }
  longacq[longacq == 0] <- 0.00001 #cannot sample in proportion of preferences if 
  # the  values are 0, the sampler doesn't work. So fix these to something very small
  for (i in 1:(b.size-1)){
    notchosen <- which(rowSums(apply(S[,c(1,2)], 1, function(a) apply(info[,c(1,2)], 1, function(b) !all(a==b)))) == b.size)
    j = i+1
    iseq <- which(rowSums(t(sapply(j:b.size, function(x) S[i,c(1,2)] == S[x,c(1,2)]))) == 2) + i
    newsize <- length(iseq)
    if(newsize !=0){
      S[iseq,c(1,2)] <- c(apply(info[sample(notchosen, newsize, prob = (longacq/sum(longacq))[notchosen]),c(1,2)], 1, as.numeric))
    }}
  S[,3] <- sapply(1:b.size, function(x) acquis[S[x,1], S[x,2]])
  return(list("NewBatch" = S[,c(1,2)]))
}

#1-D Generalized Slice Sampler without Inverse Distance Penalization 

#acquis: matrix of size J x K evaluated acquisitions (through EI or UCB)
# par.values: configuration index, (i.e., 1:100)
# month.past: past month's schedule and binary outcomes  
# b.retain: how many configurations do you want to keep from the previous month's batch
# b.size: batch size  (we set this to 8)
# toppast: is the previous month m-1's configuration INDEX 
# with the maximum posterior mean in latent preference

# More detailed annotations for this function are included in gSlice1DPenalizDist
# ignoring the sequential penalization 
# This function, instead of penalizing, simply removes previously chosen batch items from
# the candidate set 
gSlice1D <- function(acquis, par.values, month.past, b.retain = 2, b.size = 8, toppast){
  
  b.grab <- b.size - b.retain; 
  acquis[acquis < 0] <- log(1 + exp(acquis[acquis < 0]))
  J = length(par.values)
  info <- cbind(par.values, acquis); uniquepast <- unique(month.past[,c(3)])
  #slice sample
  if(!is.null(toppast)){ #in the case of using EI surface 
    pastacqui <- max(acquis[uniquepast][uniquepast != toppast])
    S <- matrix(NA, nrow = b.size, ncol= 2); S[1,] <- c(toppast, acquis[toppast])}
  if(is.null(toppast)){ #in the case of using UCB acquisition surface
    toppast <- which.max(acquis[uniquepast]); pastacqui <- max(acquis[uniquepast][uniquepast != toppast])
    S <- matrix(NA, nrow = b.size, ncol= 2); S[1,] <- c(toppast, acquis[toppast])
  }
  exclude <- c()
  minacqui = 0
  for(i in 2:b.size){
    if(i <= b.retain & i > 1){ # b.retain # element should be determined from only previous month info
      u1 <- runif(1, 0, pastacqui) #prefix the bottom to be 0; can't sample the last 
      #acquisition value evaluated at 0-valued previous best setting option so
      # just sample from maximum of the previous months' discredited values 
      accept = 0; while(accept == 0){ #mini sample from past month data
        toppast <- c(toppast, exclude)
        s0i <- info[sample(uniquepast[!(uniquepast %in% toppast)], size = 1, replace= F),1]
        s0 <- as.numeric(s0i)
        alphaprop <- acquis[s0]
        if(alphaprop > u1){S[i,] <- c(as.numeric(s0), alphaprop);
        accept = 1; exclude <- s0;}else{accept = 0; u1 <- u1*.85}
        pastacqui <- alphaprop}
    }
    if(i > b.retain){
      u1 <- runif(1, 0, pastacqui);
      accept = 0; while(accept == 0){
        si <- info[sample(par.values[!(par.values %in% toppast)], size = 1, replace= F),1]
        alphaprop <- acquis[as.numeric(si)]
        if(alphaprop > u1){S[i,] <- c(as.numeric(si), alphaprop);
        accept = 1}else{accept = 0; u1 <- u1*.85}
      }
      #print(sprintf("Collection of Item %s",i))
      pastacqui <- acquis[as.numeric(si)]
    }
  }
  for (i in 1:(b.size-1)){
    notchosen <- info[!info[,1] %in% S[,1],1]
    j = i+1
    iseq <- which(t(sapply(j:b.size, function(x) S[i,1] == S[x,1]))) + i
    newsize <- length(iseq)
    if(newsize != 0){
      S[iseq,1] <- as.numeric(info[sample(notchosen, newsize, prob = (acquis/sum(acquis))[notchosen]),1])
      S[iseq,2] <- acquis[S[iseq,1]]}
  }
  return(list("NewBatch" = S[,1]))
}


#################################### Function to choose sampling and acquisition type #############
# 2-D Batch Acquisition Strategy choice first

#acquisition types:
#"qEI" = sequential EI <-- no penalization/sampling for this
#"EI" = first step EI, used with a sampling method (or penalization)
#"UCB" = UCB
#if type = "None" and sampling = "Thompson", then Thompson sampling (sampling proportional to probability of setting being optimal)

#Sampling or Penalization methods#
#if sampling = "Penaldist" then we use distance penalized slice sampling where the acquisition function is 
# penalized by the inverse euclidean distance of the previously chosen batch element and its acquisition value 
#if sampling = "Slice" then we slice sample for each batch element

#b.retain = how many batch elements you would want to retain from the previous month (Settings)
#b.size = size of batch to go into the next month
#alpha.samps = the samples from the latent preference surface as a matrix (nrow = number of posterior samples and 
# ncol = number of settings)
# x1.values = 1:J; x2.values = 1:K; month.past = last month's schedule

acquisChoice2D <- function(type = "qEI", sampling = "None", alpha.samps, x1.values, x2.values, month.past, b.retain = 2, b.size =8){
  
  one.d.sched <- Convert2Dto1D(month.past, referenceLattice(length(x1.values),length(x2.values)))
  past <- unique(one.d.sched[,2])
  if(type == "qEI"){
    batch = rep(NA,b.size); J = ncol(alpha.samps)
    batch[1] = past[which.max(colMeans(alpha.samps[,past]))] #Restrict first two selections to previous batch!
    batch[2] = past[which.max(colMeans(alpha.samps[,past]-alpha.samps[,batch[1]]))]
    if(batch[2] == batch[1]){batch[2] = past[order(colMeans(alpha.samps[,past] - alpha.samps[,batch[1]]), decreasing = T)][2]}
    surf <- list(); batchacc <- list()
    for(b in 3:b.size){EISurface <- (colMeans(sapply(1:J,function(x) pmax(alpha.samps[,x]-apply(cbind(alpha.samps[,batch[1:(b-1)]]),1,max),0))))
    batch[b] = which.max(EISurface)
    surf[[b]] <- EISurface 
    batchacc[[b]] <- batch[1:b-1]}
    sampling <- "None"
    NewBatch <- oneDto2D(batch, x1.values, x2.values) #spits out 2-D treatment regime
  }
  if(type == "EI"){
    toppast <- past[which.max(colMeans(alpha.samps[,past]))]
    acquis <- colMeans(pmax(alpha.samps-alpha.samps[,toppast], 0))
    acquis <- matrix(acquis, ncol = length(x2.values), nrow = length(x1.values))
    if(sampling == "Penaldist"){
      NewBatch <- gSlice2DPenalizDist(acquis, x1.values, x2.values, month.past, b.retain, b.size, toppast)[[1]]
    }
    if(sampling == "Slice"){
      NewBatch <- gSlice2D(acquis, x1.values, x2.values, month.past, b.retain, b.size, toppast)[[1]]
    }
  }
  if(type == "UCB"){
    acquis <- colMeans(alpha.samps) + 1*apply(alpha.samps, 2, sd)
    acquis <- matrix(acquis, ncol = length(x2.values), nrow = length(x1.values))
    if(sampling == "Penaldist"){
      NewBatch <- gSlice2DPenalizDist(acquis, x1.values, x2.values, month.past, b.retain, b.size, toppast = NULL)[[1]]
    }
    if(sampling == "Slice"){
      NewBatch <- gSlice2D(acquis, x1.values, x2.values, month.past, b.retain, b.size, toppast = NULL)[[1]]
    }
  }
  ## If not qEI, use penalization or sampling method to obtain next 8 batches ##
  if(sampling == "Thompson" & type == "None"){
    NewBatch <- numeric(b.size)
    probofmax <- table(factor(apply(alpha.samps, 1, which.max),levels = 1:ncol(alpha.samps)))/dim(alpha.samps)[1]
    #probability of being maximum is the line above
    forsamp <- data.frame("settings" = names(probofmax), "probs"= as.numeric(probofmax)); 
    avail <- which(forsamp$probs != 0) #these are the remaining settings with non-zero probs
    #first scenario: we have too few non-zero probabilities for our previous batch elements
    forsamp$probs[forsamp$probs == 0] <- 1/(length(forsamp$settings) - length(avail))
    retained <- as.numeric(as.character(sample(forsamp[past,1], b.retain, forsamp[past,2], replace = F)))
    NewBatch[1:b.retain] <- retained #sampling from past values
    NewBatch[(b.retain+1):b.size] <- sample(forsamp[-c(retained),1], (b.size - b.retain), prob = forsamp[-c(retained),2], replace = F)
    #sampling from rest of the settings
  }
  if(type != "None"){
    NewBatch <- TwoDto1D(NewBatch, x1.values, x2.values)}
  return(list("NewBatch" = NewBatch))
}

# All same inputs as acquisChoice2D, but instead of x1.values and x2.values, 
# input the par.values, or the vector of configuration indices (i.e., 1:100 for 100 configurations)
acquisChoice1D <- function(type = "qEI", sampling = "None", alpha.samps, par.values, month.past, b.retain = 2, b.size =8){
  #qEI does sequential EI
  #If qEI, then sampling can be skipped since it will choose the top 8 from qEI; if not then specify UCB and then sampling strategy
  past <- unique(month.past[,3])
  if(type == "qEI"){
    batch = rep(NA,b.size); J = ncol(alpha.samps)
    batch[1] = past[which.max(colMeans(alpha.samps[,past]))] #Restrict first two selections to previous batch!
    batch[2] = past[which.max(colMeans(alpha.samps[,past]-alpha.samps[,batch[1]]))]
    if(batch[2] == batch[1]){batch[2] = past[order(colMeans(alpha.samps[,past] - alpha.samps[,batch[1]]), decreasing = T)][2]}
    surf <- list(); batchacc <- list()
    for(b in 3:b.size){EISurface <- (colMeans(sapply(1:J,function(x) pmax(alpha.samps[,x]-apply(cbind(alpha.samps[,batch[1:(b-1)]]),1,max),0))))
    batch[b] = which.max(EISurface)
    surf[[b]] <- EISurface 
    batchacc[[b]] <- batch[1:b-1]}
    sampling <- "None"
    NewBatch <- batch
  }
  if(type == "EI"){
    toppast <- past[which.max(colMeans(alpha.samps[,past]))]
    acquis <- colMeans(pmax(alpha.samps-alpha.samps[,toppast], 0))
    if(sampling == "Penaldist"){
      NewBatch <- gSlice1DPenalizDist(acquis, par.values, month.past, b.retain, b.size, toppast)[[1]]
    }
    if(sampling == "Slice"){
      NewBatch <- gSlice1D(acquis, par.values, month.past, b.retain, b.size, toppast)[[1]]
    }
  }
  if(type == "UCB"){
    acquis <- colMeans(alpha.samps) + 1*apply(alpha.samps, 2, sd)
    if(sampling == "Penaldist"){
      NewBatch <- gSlice1DPenalizDist(acquis, par.values, month.past, b.retain, b.size, toppast = NULL)[[1]][,1]
    }  
    if(sampling == "Slice"){
      NewBatch <- gSlice1D(acquis, par.values, month.past, b.retain, b.size, toppast = NULL)[[1]]
    }    
  }
  if(sampling == "Thompson" & type == "None"){
    NewBatch <- numeric(b.size)
    probofmax <- table(factor(apply(alpha.samps, 1, which.max),levels = 1:ncol(alpha.samps)))/dim(alpha.samps)[1]
    #probability of being maximum is the line above
    forsamp <- data.frame("settings" = names(probofmax), "probs"= as.numeric(probofmax)); 
    avail <- which(forsamp$probs != 0) #these are the remaining settings with non-zero probs
    #first scenario: we have too few non-zero probabilities for our previous batch elements
    forsamp$probs[forsamp$probs == 0] <- 1/(length(forsamp$settings) - length(avail))
    retained <- as.numeric(as.character(sample(forsamp[past,1], b.retain, forsamp[past,2], replace = F)))
    NewBatch[1:b.retain] <- retained #sampling from past values
    NewBatch[(b.retain+1):b.size] <- sample(forsamp[-c(retained),1], (b.size - b.retain), prob = forsamp[-c(retained),2], replace = F)
    #sampling from rest of the settings
  }
  return(list("NewBatch" = NewBatch))
}

