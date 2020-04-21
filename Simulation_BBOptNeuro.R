# 1-D and 2-D Simulations for Batch Bayesian Optimization Design for Optimizing a Neurostimulator
# By Adam Kaplan and Thomas Murray (2019)
# University of Minnesota, Twin Cities 

# If device calibration over a 1-D parameter space is desired then leave code as is
# If device calibration over a 2-D parameter space is desired then 
# comment out all 1-D related calls (truealphas, onesim, and par.values)
# and then uncomment the lines with 2-D function calls and objects

# WARNING: One 2-D year-long calibration simulation takes about an hour, if J = 10 and K=11, 
# we HIGHLY recommend the code in the for-loop to be parallelized 
# we used the function "foreach" from the parallel/foreach packages in R

source("functions_BBOptNeuro.R")

# 1-D simulation -- The 1-D simulations do not take too long to complete#
par.values <- 1:100
true.alpha <- pref.funcval(par.values, cores = 1, mus = c(.4), vars = c(.03), p = 1) # simulate 1-D true preference surface
# for 100 configurations
# specify a desired convergence interval 
mu = 0.5

# 1-D surface used in the article
png("TrueAlphaPriorUnimod.png")
plot(par.values,true.alpha, ylim = c(-1,1.5), type = 'b', pch = 20, main = "True Latent Parameter Surface", xlab = "Parameter Setting Values", ylab = "True Latent Preference Value")
dev.off()

# One Simulation of a Year Long Trial using a batch size of 8 with qEI as the batch acquisition strategy #
set.seed(1000)
onetrial <- onesim(par.values, true.alpha, b.size=8, b.retain=2,
                   days = 30, Qtype = 1, lambda = 0.00001,
                   type = "qEI", sampling = "None")

#2-D Simulation, IMPORTANT: comment out all 1-D parameters specified, and the 1-D sim code

#J = 10; K=11; x1.values = 1:J; x2.values = 1:K

#true.alpha <- sample2DGaussUni(J, K, c(2,2), matrix(c(2.5, -.5, -.5, 2.5), nrow = 2, ncol = 2, byrow = T))
#true.alpha <- true.alpha - mean(true.alpha)

#png("2DUnimodalNeutrality.png", width = 650, height = 200)
#  contour(x=1:K, y=1:J, z=truealpha,xlab="Frequency",ylab="Pulsewidth",
#          zlim = c(-2,2), cex.lab = 1.5, cex.axis = 1.3, las = 1,
#          nlevels = 25)
#dev.off()

#onetrial <- onesim2D(x1.values = x1.values, x2.values = x2.values, true.alpha, b.size = 8,
#b.retain = 2, days = 30, model.type = 1, type = "qEI", sampling = "None")


# Matrix of simulation metrics used in the article 
#column 1: estimated best configuration 
# column 2: at which month 
# column 3: RMSE 
# column 4: Optimal Preference Difference
# column 5: Trial Quality 
onetrial$SimMeasures # matrix of length 12 by 5, 12 being the number of months, 5 for metrics

#Type 1 and Type 2 Brier Scores 
# Brier scores closer to 0 are optimal 
onetrial$BrierTraject

#Early stopping probabilities for preference neutrality and calibration convergence 
#column vector of length of 12 (for 12 months)
onetrial$FutilityMat # Preference Neutrality probabilities 
onetrial$SuperiorMat # Calibration Convergence probabilities 

# These results indicate that at no time during the trial we would stop the trial for preference neutrality with 
# a cut off of p_0 = 15%
# however, we may stop the trial for calibration convergence at around month 2 if we specified a tolerance 
# s_0 = 80%. If you require "more evidence" then increase s_0 to stop the trial for calibration convergence

onetrial$SimMeasures[min(which(onetrial$SuperiorMat > .8)),]
# This obtains the proper calibration measures for the month we stopped the calibration if we 
# wanted to stop for calibration convergence with a s_0 = 80%. 

###########################################################################
###########################################################################
########## Conducting 100 Simulations of 1-D (or 2-D) Configuration #######
########## for q-EI batch acquisition function only #######################
###########################################################################
###########################################################################

# Grid of Preference Neutrality and Calibration Convergence values
p.cuts <- seq(0, .2, by = .005); s.cuts <- seq(.6, 1, by= .005)

times <- 100

output <- list()
for(i in 1:times){
  
  set.seed(i)
  onetrial <- onesim(par.values = 1:100, true.alpha = true.alpha,
                     b.size = 8, b.retain = 2, days = 30, Qtype = 1, lambda = 0.0001,
                     type = "qEI", sampling = "None")
  
  # If 2-D is desired, comment out the onetrial line above #
  #onetrial <- onesim2D(x1.values = x1.values, x2.values = x2.values, true.alpha, b.size = 8,
  #b.retain = 2, days = 30, model.type = 1, type = "qEI", sampling = "None")
  
  futilsmats <- t(sapply(1:length(p.cuts), function(s)onetrial$FutilityMat > p.cuts[s]))
  futilsmats <- apply(futilsmats, 1, function(l){o <- min(which(l));
  o[o == Inf] <- 12; return(o)})
  
  supersmats<- t(sapply(1:length(s.cuts), function(s)onetrial$SuperiorMat > s.cuts[s]))
  supersmats <- apply(supersmats, 1, function(l){o <- min(which(l));
  o[o == Inf] <- 12; return(o)}) 
  
  futDat <- t(sapply(1:length(futilsmats), function(x) onetrial[[1]][futilsmats[x],]))
  #matrix of row-length length(p.cuts) and column-length of number of metrics 
  supDat <- t(sapply(1:length(supersmats), function(x) onetrial[[1]][supersmats[x],]))
  #matrix of row-length length(s.cuts) and column-length of number of metrics 
  output[[i]] <- list("FutilityMat" = futDat, "SuperiorityMat" = supDat)
  #print(sprintf("Just completed simulation of trial %d", i))
  
}

################# Doing this in parallel is a lot faster #####################

################## Reporting results without stopping the trial early for q-EI ################################
ress <- t(sapply(1:length(output), function(x) output[[x]]$SuperiorityMat[nrow(output[[x]]$SuperiorityMat),]))

numericres <- rbind(apply(ress[,c(5,4,3)], 2, mean),
                                       apply(ress[,c(5,4,3)], 2, function(y) sqrt(var(y)/nrow(ress))))

# Trial Quality, Optimal Preference Difference, and RMSE 
numericres

#Only if 1-D is present, this would make sense
plot(density(ress[,1]), main = "Density of Selected Configurations", xlab = "Configuration Index", las = 1)

############################## Now obtaining results for early stopping################################### 
p.cut4me <- which(p.cuts == .15)
s.cut4me <- which(s.cuts == .80)

###################Only early stopping for preference neutrality with a cut off of .15
ress1 <- t(sapply(1:length(output), function(x) output[[x]]$FutilityMat[p.cut4me,]))

numericres <- cbind(rbind(apply(ress1[,c(5,4,3)], 2, mean),
                          apply(ress1[,c(5,4,3)], 2, function(y) sqrt(var(y)/nrow(ress1)))),
                    c(mean(ress1[,c(2)]), sd(ress1[,c(2)])))

numericres
###################Only early stopping for Calibration Convergence with a cut off of .80
ress2 <- t(sapply(1:length(output), function(x) output[[x]]$SuperiorityMat[s.cut4me,]))

numericres <- cbind(rbind(apply(ress2[,c(5,4,3)], 2, mean),
                          apply(ress2[,c(5,4,3)], 2, function(y) sqrt(var(y)/nrow(ress2)))),
                    c(mean(ress2[,c(2)]), sd(ress2[,c(2)])))

numericres
####################For both: see which stopping rule - month is earlier for each simulation ##################

indics <- ress1[,2] < ress2[,2];# specifying the month of stopping that is earlier between the two rules 
ress3 <- rbind(ress1[indics,], ress2[!indics,]) # putting them altogether for summary measures 

numericres <- cbind(rbind(apply(ress3[,c(5,4,3)], 2, mean),
                    apply(ress3[,c(5,4,3)], 2, function(y) sqrt(var(y)/nrow(ress3)))),
                    c(mean(ress3[,c(2)]), sd(ress3[,c(2)])))

# Trial Quality, Optimal Preference Difference, and RMSE, and mean(sd) month of stopping
numericres


