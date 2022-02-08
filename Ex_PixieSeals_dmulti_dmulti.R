#************************#
# EXAMPLE 1: PIXIE SEALS #
#************************#

# Population dynamics: 
# Stage 1 = age 1 subadult --> conditional on survival, may mature (-> 4) or not (-> 2)
# Stage 2 = age 2 subadult --> conditional on survival, may mature (-> 4) or not (-> 3)
# Stage 3 = age 3 subadult --> matures conditional on survival (-> 4)
# Stage 4 = mature adult --> remains in stage conditional on survival, reproduces

# (Stage 5 = dead --> absorbing state)


# Data collection:
# The whole population is counted (with Poisson-distributed error) at each time-step
# Data does note allow distinguishing between stages

mySeed <- 88
set.seed(mySeed)

# Simulation parameters #
#-----------------------#

## Survival rates
S <- c(0.7, 0.8, 0.85, 0.95)

## Maturation rates
pMat <- c(0.1, 0.5)

## Reproductive rate
R <- 0.4

## Initial population size
Ntot.init <- 1000

## Initial population structure
SS.init <- c(0.3, 0.2, 0.1, 0.4) 

## Number of stages 
Amax <- 4

## Number of simulation steps (years)
Tmax <- 5


# Population dynamics simulation #
#--------------------------------#

## Assemble transition probability matrix (time-invariant)
Psi <- matrix(NA, nrow = Amax, ncol = Amax+1)
Psi[1,] <- c(0, S[1]*(1-pMat[1]), 0, S[1]*pMat[1], 1-S[1])
Psi[2,] <- c(0, 0, S[2]*(1-pMat[2]), S[2]*pMat[2], 1-S[2])
Psi[3,] <- c(0, 0, 0, S[3], 1-S[3])
Psi[4,] <- c(0, 0, 0, S[4], 1-S[4])

## Matrix for storing population size/structure
N <- matrix(NA, nrow = Amax, ncol = Tmax+1)

## Array for storing population transitions
transN <- array(NA, dim = c(Amax, Amax+1, Tmax))

## Simulate initial population
N[,1] <- rmultinom(1, size = Ntot.init, prob = SS.init)

## Simulate population dynamics over time
for(t in 1:Tmax){
  
  # Reproduction
  N[1,t+1] <- rbinom(1, size = N[4,t], prob = R)
  
  # Survival & stage transitions
  for(a in 1:Amax){
    transN[a, 1:(Amax+1), t] <- rmultinom(n = 1, size = sum(N[a,t]), prob = Psi[a, 1:(Amax+1)])
  }
  
  # Summing to population size
  N[2:Amax, t+1] <- colSums(transN[1:Amax, 2:Amax, t])
  
}


# Data collection simulation #
#----------------------------#

## Vector to store observations
obs.Ntot <- rep(NA, Tmax+1)

## Simulate data collection
for(t in 1:(Tmax+1)){
  obs.Ntot[t] <- rpois(1, lambda = sum(N[,t]))
}


# Formulate Nimble model for estimating population sizes #
#--------------------------------------------------------#

## Assumptions: 
#   - True population size is unknown = latent quantity of interest
#   - Demographic rates and initial stage proportions are known without error


## Model code
seal.code <- nimbleCode({
  
  # Population dynamics model #
  #---------------------------#
  
  # Initializing population model
  N[1:Amax,1] ~ dmulti(size = Ntot[1], prob = SS.init[1:Amax])
  
  # Population dynamics
  for(t in 1:Tmax){
    
    # Reproduction
    N[1, t+1] ~ dbin(size = N[4,t], prob = R)
    
    # Survival & stage transitions
    for(a in 1:Amax){
      transN[a, 1:(Amax+1), t] ~ dmulti(size = N[a, t], prob = Psi[a, 1:(Amax+1)])
    }
    
    # Summing to population size
    for(a in 2:Amax){
      N[a, t+1] <- sum(transN[1:Amax, a, t])
    }
  }
  
  
  # Data likelihood(s) #
  #--------------------#
  
  # Count data likelihood
  for(t in 1:(Tmax+1)){
    obs.Ntot[t] ~ dpois(Ntot[t])
  }
  
  
  # Priors and constraints #
  #------------------------#
  
  Ntot[1] <- round(cont.Ntot1)
  cont.Ntot1 ~ dunif(0, 2000)
  
  for(t in 2:(Tmax+1)){
    Ntot[t] <- sum(N[1:Amax,t])
  }
})


# Estimate population sizes using Nimble #
#----------------------------------------#

## Arrange data and constants
test.data <- list(obs.Ntot = obs.Ntot, SS.init = SS.init)
test.constants <- list(Amax = Amax, Tmax = Tmax,
                       R = R, Psi = Psi)

## Set Initial values 
init.fun <- function(){
  list(cont.Ntot1 = rpois(1, Ntot.init))
}
  
Inits <- init.fun()

## Set parameters to monitor
params <- c('N', 'Ntot')

## Test run
test.run <- nimbleMCMC(code = seal.code, 
                       constants = test.constants, 
                       data = test.data, 
                       inits = Inits, 
                       monitors = params, 
                       niter = 100, nburnin = 0, nchains = 1, thin = 1, 
                       setSeed = 0, samplesAsCodaMCMC = TRUE)

test.run

# This illustrates the original problem: 
# The sampler fails to update latent multinomial nodes (with the exception of
# nodes for the terminal time-step, here t=6)


