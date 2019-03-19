# "A Preferential Attachment Model for the Stellar Initial Mass Function" 
#  by Cisewski-Kehe, Weller, and Schafer

# Generate masses of stars in stellar cluster via Preferential Attachment mechanism

clusterSim <- function(Mtot, theta){	
	#Mtot = total mass of cluster to be converted to stars
	#theta = (k0, alpha, gamma)
	#k0 = 1/lambda, controls size of mass components (k0 > 0)
	#alpha = probability of a new star forming (0 < alpha < 1)
	#gamma = growth component (gamma>0, 1 is linear attachment)
	k0 <- theta[1] 
	alpha <- theta[2]
	gamma <- theta[3]
	
#Initiate cluster
Mobs <- rexp(1,1/k0)

#Cluster growth
while(sum(Mobs) < Mtot){
	if(runif(1) > alpha){
		pVec <- Mobs^gamma
		Mobs <- Mobs + rexp(1,1/k0) * as.vector(rmultinom(1,1,pVec)) 
	}else{
		Mobs <- c(Mobs, rexp(1,1/k0)) #new star
		}
	}

# Output
out <- list()
out$Mimf <- Mobs
return(out)
}



