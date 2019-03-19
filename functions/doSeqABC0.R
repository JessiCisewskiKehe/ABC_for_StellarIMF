# "A Preferential Attachment Model for the Stellar Initial Mass Function" 
#  by Cisewski-Kehe, Weller, and Schafer

# Iterative steps of the ABC algorithm

doSeqABC <- function(nIter = 1, input, tol, which.tol, inputDist, prev.draws, prev.weights){
	# See doInitABC1.R for explanations of inputs
Mtot <- input$Mtot	
Mobs <- input$M
age <- input$age
sigma <- input$sigma
obsProb <- input$obsProb
thetaTrue <- input$thetaTrue
which.input <- input$which.input
parameter.grid <- input$parameter.grid

i <- 0
draws <- 0
while(i < nIter){
	draws <- draws + 1
	prev.wmean <- t(prev.draws)%*%prev.weights
	prev.center <- prev.draws - matrix(rep(prev.wmean,Ndraws),byrow = T, nrow = Ndraws)
	var.prev <- t(prev.center)%*%diag(prev.weights)%*%prev.center

# Draw from the previous values
	index.samp <- sample(1:Ndraws, 1, prob = prev.weights)
	theta.gen <- if(length(which.input) >1){
		prev.draws[index.samp,]
		}else prev.draws[index.samp]


# Perturb it
	thetaProp <- thetaTrue
	thetaProp[which.input] <- rmvnorm(1, theta.gen, 2*var.prev)
# Try again if it's outside the prior ranges
	T1 <- thetaProp[1] < parameter.grid[1,1] | thetaProp[1] > parameter.grid[2,1] 
	T2 <- thetaProp[2] < parameter.grid[1,2] | thetaProp[2] > parameter.grid[2,2]
	T3 <- thetaProp[3] < parameter.grid[1,3] | thetaProp[3] > parameter.grid[2,3]
	while(T1 | T2 | T3){
		thetaProp[which.input] <- rmvnorm(1, theta.gen, 2*var.prev)
		T1 <- thetaProp[1] < parameter.grid[1,1] | thetaProp[1] > parameter.grid[2,1] 
		T2 <- thetaProp[2] < parameter.grid[1,2] | thetaProp[2] > parameter.grid[2,2]
		T3 <- thetaProp[3] < parameter.grid[1,3] | thetaProp[3] > parameter.grid[2,3]
	}

	# Generate simulated stars
	clusterOut = clusterSim(Mtot, theta = thetaProp)
	Mimf <- clusterOut$Mimf

	#Put through forward model
	outForward = forwardModel(Mimf, age, sigma, obsProp)
	Mprop <- outForward$Mprop

	Nprop <- length(Mprop)
	MmaxProp <- max(Mprop)
	
	dist = if(length(Mprop) >3){dist = distance.function(Mprop, Mobs, inputDist)
		}else{dist= tol+1}
	
	if(sum(dist[which.tol]<=tol[which.tol])==length(which.tol)){i <- i + 1}


}
# Get the density value of the chosen draw
theta.kern <- 1/sum(prev.weights*apply(matrix(prev.draws,nrow = Ndraws), 1, dmvnorm, x=thetaProp[which.input], sigma = 2*var.prev, log=FALSE))

out <- list()
out$draws <- draws
out$N <- Nprop
out$Mmax <- MmaxProp
out$theta <- thetaProp[which.input]
out$kernel <- theta.kern
out$D <- dist
out$Mprop <- Mprop
return(out)
}
