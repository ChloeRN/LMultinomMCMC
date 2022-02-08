library(nimble)
library(gtools)

set.seed <- 0


# 0) Set constants
#-----------------

Zo <- 3 # Number of origin states
#Zt <- Zo*3 # Number of target states
Zt <- Zo+1 # Number of target states
Tmax <- 3 # Number of time-steps


# 1) Simulate data
#-----------------

# 1.1) Transition probabilities
U <- array(0, dim = c(Zo, Zt, Tmax))

for(t in 1:Tmax){
	
	for(z in 1:Zo){
		U[z, 1:Zt, t] <- rdirichlet(1, rep(1, Zt)) 
	}
}


# 1.2) Starting numbers
Ninit <- round(matrix(runif(Zo*Zt, 100, 1000), nrow = Zo, ncol = Zt))

# 1.3) Stochastic realizations
N <- array(NA, dim = c(Zo, Zt, Tmax+1))

N[,,1] <- Ninit

for(t in 1:Tmax){
	
	for(z in 1:Zo){
		N[z, 1:Zt, t+1] <- rmultinom(n = 1, size = sum(N[1:Zo, z ,t]), prob = U[z, 1:Zt, t])
	}
}


# 2) NIMBLE model code
#---------------------

dummy.code <- nimbleCode({
  
  
  ## Process model (multinomial sampling)
  
  for(t in 1:Tmax){
    

    #N[1, 1:Zt, t+1] ~ dmulti(p = U[1, 1:Zt, t], size = sum(N[1:Zo, 1,t]))    
    
    for(z in 1:Zo){
    	N[z, 1:Zt, t+1] ~ dmulti(p = U[z, 1:Zt, t], size = sum(N[1:Zo, z, t]))
    }
    
  } 
  
  ## Priors for initial numbers

  for(z in 1:Zo){
  	for(y in 1:Zt){
  		initN[z,y] ~ dunif(0, 1000)
  		N[z,y,1] <- round(initN[z,y])
  	}  
  }
  
})


# 3) Test run
#------------

# 3.1) Arrange data and constants
dummy.data <- list(U = U)
dummy.constants <- list(Zo = Zo, Zt = Zt, Tmax = Tmax)

# 3.2) Initial values 

# Starting numbers only (preferred in practice)
Inits <- list(initN = round(matrix(runif(Zo*Zt, 100, 1000), nrow = Zo, ncol = Zt)))

# All numbers (to check)
#Inits <- list(initN = Ninit, N = N)


# 3.3) Test run
test.run <- nimbleMCMC(code = dummy.code, constants = dummy.constants, data = dummy.data, inits = Inits, monitors = c('N'), niter = 100, nburnin = 0, nchains = 1, thin = 1, setSeed = 0, samplesAsCodaMCMC = TRUE)

test.run

# Run with nimble 0.9.0 on Mac: 
#saveRDS(test.run, file = 'MultinomEx_v0.9.0_Mac.rds')
# --> Nodes are fixed / not updating (original problem)

# Run with nimble 0.12.1 on Mac: 
#saveRDS(test.run, file = 'MultinomEx_v0.12.1_Mac.rds')
# --> Nodes update

# Run with nimble 0.12.1 on Windows: 
#saveRDS(test.run, file = 'MultinomEx_v0.12.1_Win.rds')
# --> Nodes update

# SIDE NOTE: I was unable to test with nimble 0.9.0 on Windows due to some
#            library issues. 
