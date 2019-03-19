# "A Preferential Attachment Model for the Stellar Initial Mass Function" 
#  by Cisewski-Kehe, Weller, and Schafer

# Initial step of the ABC algorithm

doInitABC <- function(nIter = 1, input, tol, which.tol, inputDist){
# INPUTS
#	nIter = number of iterations (ABC draws) desired ("1" is the only value that makes sense)
#	input = object storing Mtot, M, age, sigma, obsProb, thetaTrue (when available in simulation settings), which.input, parameter.grid
#	tol = distance threshold
#   which.tol = specify the index of the summary statistics/distances used
#   inputDist = used in distance function (usually KDE of observed IMF)

Mtot <- input$Mtot	
Mobs <- input$M
age <- input$age
sigma <- input$sigma
obsProb <- input$obsProb
thetaTrue <- input$thetaTrue
which.input <- input$which.input
parameter.grid <- input$parameter.grid

i <- 0
#########################################################################	
	# Prior distributions
#########################################################################	
	thetaProp <- thetaTrue
	thetaProp[which.input] <- runif(length(which.input),parameter.grid[1,which.input], parameter.grid[2, which.input])

# Generate simulated stars
clusterOut = clusterSim(Mtot, theta = thetaProp)
Mimf <- clusterOut$Mimf

#Put through forward model
outForward = forwardModel(Mimf, age, sigma, obsProb)
Mprop <- outForward$Mprop

Nprop <- length(Mprop)
MmaxProp <- max(Mprop)

if(length(Mprop)>3){
dist = distance.function(Mprop, Mobs, inputDist)
} else{
	dist = rep(1^3,length(which.tol))
	}

out <- list()
out$Mtot <- Mtot
out$N <- Nprop
out$Mmax <- MmaxProp
out$theta <- thetaProp[which.input]
out$D <- dist
out$Mprop <- Mprop
return(out)
}
