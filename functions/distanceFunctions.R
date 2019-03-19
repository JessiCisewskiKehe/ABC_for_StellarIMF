# "A Preferential Attachment Model for the Stellar Initial Mass Function" 
#  by Cisewski-Kehe, Weller, and Schafer

# Distance functions for ABC algorithm



# Shape of IMF
kde.distance = function(x,y, inputDist){
	#x = star masses from simulated cluster
	#y = not used
	#inputDist = contains KDE information about observed cluster
	xEval <- inputDist$xEval
	ObsDens <- inputDist$density 

	kdeProp <- density(log10(x),bw='nrd',kernel='gaussian')
	densProp <- approx(kdeProp$x,kdeProp$y,xout=xEval)$y
	densProp[is.na(densProp)==TRUE]<-0
	
	# Trapezoidal rule to approximate the integral numerically
	traps <- (densProp - ObsDens)^2
	dist <- 0.5*sum((xEval[2] - xEval[1])*(traps[-1] + traps[-length(traps)]))
	
	return(dist)
}	


# Star counts
star.distance = function(x,y, inputDist){
	#x = star masses from simulated cluster
	#y = star masses from observed cluster
	#inputDist = not used
	max(abs(1 - length(x)/length(y)),abs(1 - length(y)/length(x)))
}

