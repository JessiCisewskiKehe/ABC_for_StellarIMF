# "A Preferential Attachment Model for the Stellar Initial Mass Function" 
#  by Cisewski-Kehe, Weller, and Schafer

# Forward model of the ABC algorithm

forwardModel = function(Mimf, age, sigma, obsProp){
	#Mimf = masses of stars generated 
	#sigma = measurement error
	#obsProp = function specifying the probability of observing a star (i.e., the completeness function)
	
Mprop <- Mimf

#Aging effects
if(age > 0){
	# Mprop <- Mprop[Mprop <= age^(-2/5)*10^(8/5)]
	Mprop <- Mprop[Mprop <= age^(-1/3)*10^(4/3)]
	}

#Completeness function
if(length(Mprop) > 1){
	pObs <- obsProb(Mprop)
	ObsVal <- apply(matrix(pObs,ncol=1),1,rbinom,n=1,size=1)
	Mprop <- Mprop[ObsVal==1]
}

# Lognormal uncertainties
if(length(Mprop) > 1 & sigma>0){
	Mprop <- exp(log(Mprop) + rnorm(length(Mprop),0,sigma))
	}

# Output
out <- list()
out$Mprop <- Mprop
return(out)
}



