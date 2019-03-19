
Code for the ABC algorithm proposed in "A Preferential Attachment Model for the Stellar Initial Mass Function" by Cisewski-Kehe, Weller, and Schafer

# /functions

clusterSim0.R:  generate cluster masses using preferential attachment

distanceFunctions.R: distance functions for ABC algorithm

doInitABC1.R:  initial step of ABC algorithm

doSeqABC0.R: iterative steps of ABC algorithm

forwardModel0.R:  forward model of ABC algorithm



# Required R packages:
library(parallel)  #code is set to run in parallel
library(mvtnorm)
library(ks)


# example1.R:  uses data simulated from the proposed generative model with specified parameter values.  A linear ramp completeness function is used, the cluster is aged 30 Myr, and lognormal measurement error is applied.  User inputs can be updated.

# example2.R:  uses kroupa_model.txt masses simulated from the Kroupa (2001).  No observational effects (i.e. no completeness, aging, or measurement error) have been applied.