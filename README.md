# LMultinomMCMC
Repository for working on Nimble samplers for models with latent multinomial likelihoods.


## Pixie Seal Example
Pixie seals (concept art may follow later) are ethereal pinnipeds with a life history roughly inspired by ringed seals. 
The example is a model for describing the stage-structured population dynamics of pixie seals based on an annual population census. 

At census, all individuals within the population are counted but distinguishing life stages is impossible by naked eye due to the pixie seals' natural luminescence. They also like to gather in groups, resulting in some counting error. We assume this error can be appropriately represented by a Poisson distribution, and hence model the likelihood of the annual count data (`obs.Ntot[t]`) as 
```
  obs.Ntot[t] ~ dpois(Ntot[t])
```
where `Ntot[t]` is the latent population size in year ´t´, summed over all age classes, i.e. `Ntot[t] <- sum(N[1:Amax, t])` (`Amax` = number of stages). 

While life stages cannot be distinguished at census, they matter a lot for population dynamics since survival (`S[a]`) and maturation rates (`pMat[a]`) are age-dependent, and only mature individuals produce any offspring. The population model is therefore structured by four distinct stages: 

`a=1`: Age 1 subadults, cannot reproduce, may mature with probability `pMat[1]` \
`a=2`: Age 2 subadults, cannot reproduce, may mature with probability `pMat[2]` \
`a=3`: Age 3 subadults, cannot reproduce, will mature (probability = 1) \
`a=4`: Mature adults, can reproduce

The first multinomial likelihood in the pixie seal model appears when determining the stage structure of the inital population size (t=1): 
```
  N[1:Amax, 1] ~ dmultinom(size = Ntot[1], prob = SS.init[1:Amax])
```
The proportions of the population in each stage (`SS.init[a]`) are assumed to be known a priori for the sake of the example. 

The second multinomial likelihood forms the core of the population projection, describing the movement of pixies seals among life stages over time:
```
  # Survival & stage transitions
  for(a in 1:Amax){
    transN[a, 1:(Amax+1), t] ~ dmulti(size = N[a, t], prob = Psi[a, 1:(Amax+1)])
  }
    
  # Summing to population size
  for(a in 2:Amax){
    N[a, t+1] <- sum(transN[1:Amax, a, t])
  }
```
Here, we the array `transN` stores information on how many individuals transition from any stage to any other stage over the time interval t to t+1. The transition matrix `Psi` contains transition rates consisting of stage-specific survival (`S[a]`) and maturation rates (`pMat[a]`); see code for details. \
The number of stage 1 individuals in the next year is the outcome of reproduction of mature individuals:
```
  N[1, t+1] ~ dbin(size = N[4,t], prob = R)
```

The multinomial likelihoods in the pixie seal example illustrate the current issue we have with Nimble's samplers for nodes which have a latent "size" argument. 
The initial four code files illustrate this by means of implementing the pixie seal model with both, only one, or none of the multinomial likelihoods:
- `Ex_PixieSeals_dmulti_dmulti.R`: model using multinomial likelihood for initializing and projecting population model
- `Ex_PixieSeals_dbin_dmulti.R`: model using multinomial likelihood only for projecting population model
- `Ex_PixieSeals_dmulti_dbin.R`: model using multinomial likelihood only for initializing population model
- `Ex_PixieSeals_dbin_dbin.R`: model without multinomial likelihoods



