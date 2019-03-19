# "A Preferential Attachment Model for the Stellar Initial Mass Function" 
#  by Cisewski-Kehe, Weller, and Schafer

# Example using "observations" simulated from forward model

rm(list=ls(all=TRUE))  # clear all variables
library(parallel)  #code is set to run in parallel
library(mvtnorm)
library(ks)

###--------------------USER INPUT NEEDED
theta.use <- c(0.25, .3, 1) #specify desired values for (1/lambda, alpha, gamma)

PathToFunctions <- getwd() #You may need to update the directory
PathToFunctions <- paste(PathToFunctions, "/functions/", sep = "")

nCores <- detectCores()
nCores

#Specify completeness function (currently set as linear ramp between C_low and C_up)
C_low <- 0  
C_up <- 1
obsProb <- Vectorize(function(x){
	if(x < C_low){return(0)}
	else if(x > C_up){return(1)}
	else{return((x-C_low)/(C_up - C_low))}
}, 'x')


#algorithm parameters
Ndraws <- 50	 #number of particles (should be larger, but requires more computational time)
Ninit <- Ndraws*10  #number of initial particles
T_Iter <- 3  #number of iterations (should be larger, but requires more computational time)
which.quantile = .25  #quantile for decreasing tolerance

#parameter values
which.input <- c(1,2,3)  # use all 3 parameters, but can fix parameter values by leaving out of this list
k0input <- theta.use[1]
alphaInput <- theta.use[2]
gammaInput <- theta.use[3]
MtotInput <- 100 #(can be larger, but requires more computational time) 
sigma <-  .25  #measurement error
age <-  30  #age of cluster in Myr
	

which.tol = c(1,2)  # using both summary statistics
tolSeq <- c(NA,NA)	
thetaTrue <- cbind(k0input, alphaInput, gammaInput) #check on alg performance
grid.k <- c(.01, 2)  #prior range
grid.alpha <- c(0, 1)  #prior range
grid.gamma <- c(0, 2)  #prior range
parameter.grid <- cbind(grid.k, grid.alpha, grid.gamma)
parameter.grid.norm <- 1/(parameter.grid[2,]-parameter.grid[1,])


#Read in functions
source(paste(PathToFunctions,"clusterSim0.R", sep = ""))  #Generates initial star masses
source(paste(PathToFunctions,"distanceFunctions.R", sep = ""))  #Distance functions
source(paste(PathToFunctions,"doInitABC1.R", sep = ""))  # Adaptive start
source(paste(PathToFunctions,"doSeqABC0.R", sep = ""))  # Sequential ABC
source(paste(PathToFunctions,"forwardModel0.R", sep = ""))  # Forward model to get MF


distance.function = function(x,y, inputDist){
	out = cbind(kde.distance(x,y, inputDist), star.distance(x,y, inputDist))
	return(out)
}




#Generate the 'observations'
sim <- clusterSim(Mtot = MtotInput, theta = thetaTrue)
Mimf <- sim$Mimf

#Put through forward model
outForward <- forwardModel(Mimf, age, sigma, obsProb)
MF <- outForward$Mprop  # these are the "real" data

# Evaluate the density 
kdeObs <- density(log10(MF),bw='nrd',kernel='gaussian')
# Evaluate the density over a grid
xEval <- seq(log10(10^-7), log10(max(MF))*1.5, 0.001)
densVals <- approx(kdeObs$x,kdeObs$y,xout=xEval)$y
densVals[is.na(densVals)==TRUE]<-0

inputDist <- list() #inputs for distance function
inputDist$xEval <- xEval
inputDist$density <- densVals

#Setup algorithm inputs
input.list <- list(M = MF, Mtot = MtotInput, Ndraws = Ndraws, tol = tolSeq, T_Iter = T_Iter, age = age, sigma = sigma, which.input = which.input, thetaTrue = thetaTrue, parameter.grid = parameter.grid, obsProb = obsProb)


# Output that we wish to store - could add more
totDraws <- NULL		# total number of prior draws
post.theta <-NULL
post.weights <- NULL
post.Mprop <- NULL
post.Mmax <- NULL



###  Adaptive Start
t_Iter <- 1
ABCdrawsInit <- mclapply(rep(1,Ninit), doInitABC, input = input.list, tol = tolSeq, which.tol = which.tol, inputDist = inputDist, mc.cores =nCores) 
cat('Initial ABC done','\n')


#Select N particles with smallest distance
Dvals <- matrix(unlist(lapply(ABCdrawsInit,function(x){return(x$D)})),ncol = length(which.tol), byrow = TRUE)
D0 <- rank(apply(scale(Dvals, center = FALSE)^2,1,sum), ties.method = "random")
pickNdraws <- D0<=Ndraws
Dvals2 <- as.matrix(Dvals[pickNdraws,])
	
	
#Pull everything together 	
theta.Current <- matrix(unlist(lapply(ABCdrawsInit,function(x){return(x$theta)})),nrow= Ninit,byrow=T)[pickNdraws, ]

weights <- rep(1/Ndraws,Ndraws)
post.weights[[t_Iter]] <- weights
post.ess <- Ndraws
post.theta[[t_Iter]] <- theta.Current
Nprop <- unlist(lapply(ABCdrawsInit,function(x){return(x$N)}))
post.Nprop <- Nprop [pickNdraws]
totDraws <- Ninit

post.tolSeq <- apply(Dvals2,2,max)
tolSeq <- apply(Dvals2, 2, quantile, which.quantile)
cat('Original tolerance', post.tolSeq,'\n')
cat('New tolerance is', tolSeq,'\n')	


# NOW ITERATE ON T
for(t_Iter in 2:T_Iter){
	print(t_Iter)
	ABCdrawsT <- mclapply(rep(1,Ndraws), doSeqABC, input = input.list, tol = tolSeq, which.tol = which.tol, inputDist = inputDist, prev.draws = theta.Current, prev.weights = weights, mc.cores = nCores) 
	theta.Current <- matrix(unlist(lapply(ABCdrawsT,function(x){return(x$theta)})),nrow=Ndraws,byrow=T)
	post.theta[[t_Iter]] <- theta.Current
	Dvals <- matrix(unlist(lapply(ABCdrawsT,function(x){return(x$D)})),ncol = length(tolSeq), byrow = TRUE)
	totDraws <- c(totDraws,sum(unlist(lapply(ABCdrawsT,function(x){return(x$draws)}))))
	Nprop <- unlist(lapply(ABCdrawsT,function(x){return(x$N)}))
	Ncomp <- unlist(lapply(ABCdrawsT,function(x){return(x$Ncomp)}))
	kern.wt <- unlist(lapply(ABCdrawsT,function(x){return(x$kernel)}))
	weights <- kern.wt/sum(kern.wt) 
	post.weights[[t_Iter]] <- weights
	ess <- 1/sum(weights^2)
	
	post.tolSeq <- rbind(post.tolSeq, tolSeq)
	tolSeq <- apply(Dvals, 2, quantile, which.quantile)

	cat('Effective sample size is',ess,'\n')
	cat('New tolerance is', tolSeq,'\n')	
}



#Plot marginal posteriors
which.step <- t_Iter
par(mfrow = c(1,3))
plot(density(post.theta[[which.step]][, 1],weights = post.weights[[which.step]]), lwd = 3, main = expression(lambda^{-1}), xlim = parameter.grid[,1])
abline(v = thetaTrue[1], col = "red",lty = 2, lwd = 3)
abline(v = parameter.grid[1,1], col = "gray",lty = 2, lwd = 3)
abline(v = parameter.grid[2,1], col = "gray",lty = 2, lwd = 3)

plot(density(post.theta[[which.step]][, 2],weights = post.weights[[which.step]]), lwd = 3, main = expression(alpha), xlim = parameter.grid[,2])
abline(v = thetaTrue[2], col = "red",lty = 2, lwd = 3)
abline(v = parameter.grid[1,2], col = "gray",lty = 2, lwd = 3)
abline(v = parameter.grid[2,2], col = "gray",lty = 2, lwd = 3)

plot(density(post.theta[[which.step]][, 3],weights = post.weights[[which.step]]), lwd = 3, main = expression(gamma), xlim = parameter.grid[,3])
abline(v = thetaTrue[3], col = "red",lty = 2, lwd = 3)
abline(v = parameter.grid[1,3], col = "gray",lty = 2, lwd = 3)
abline(v = parameter.grid[2,3], col = "gray",lty = 2, lwd = 3)



