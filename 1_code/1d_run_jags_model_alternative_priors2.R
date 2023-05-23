library(jagsUI)

load('0_data/data.RData')
str(data)


# Specify model in BUGS language
cat(file="1_code/jags_model_priors2.txt","
model{

## Likelihood

for(i in 1:nsites){        # Loop over sites
  for(k in 1:nspecies){    # Loop over species

    z[i,k] ~ dbern(psi[i,k])
    logit(psi[i,k]) <- alpha0[k]
                        + alpha1*substrate[i,k]
                        + alpha2[k]*altitude.st.2[i]
                        + alpha3[k]*altitude.st[i]
                        + alpha4[k]*precipitation.st[i]
                        + alpha5[k]*precipitation.st.2[i]
                        + alpha6[k]*solar.st[i]
                        + alpha7[k]*solar.st.2[i]

    for(j in visits[i,1]:visits[i,2]){  # Loop over 1-2 visits

      y[i,j,k] ~ dbern(z[i,k]*p[i,j,k])
      # y.sim[i,j,k] ~ dbern(z[i,k]*p[i,j,k])    # only for goodness-of-fit assessment
      logit(p[i,j,k]) <- beta1[observer[i,j],k]
                          + beta2*conspicuousness[k]
                          + beta3*identifiability[k]
                          + beta4*experience[i,j,k]

    }
  }
}

# Define random effects

for(k in 1:nspecies){

  # For the occupancy part of the model
  alpha0[k] ~ dnorm(mu.alpha0, tau.alpha0) # Random species effect, intercept
  alpha2[k] ~ dnorm(mu.alpha2, tau.alpha2) # Random species effect, slope altitude.st
  alpha3[k] ~ dnorm(mu.alpha3, tau.alpha3) # Random species effect, slope altitude.st.2
  alpha4[k] ~ dnorm(mu.alpha4, tau.alpha4) # Random species effect, slope precipitation.st
  alpha5[k] ~ dnorm(mu.alpha5, tau.alpha5) # Random species effect, slope precipitation.st.2
  alpha6[k] ~ dnorm(mu.alpha6, tau.alpha6) # Random species effect, slope solar.st
  alpha7[k] ~ dnorm(mu.alpha7, tau.alpha7) # Random species effect, slope solar.st.2

  # For the detection part of the model
  for(o in 1:nobservers){
    beta1[o,k] ~ dnorm(mu.beta1[o], tau.beta1[o]) # Random species effect on detection within each observer
  }
}


## Priors for occupancy process

mu.alpha0.prior ~ dunif(0,1)         # Mean hyperparameter on intercept
mu.alpha0 <- logit(mu.alpha0.prior)
sd.alpha0 ~ dt(0, 1/25^2, 1) T(0,) # SD hyperparameter
tau.alpha0 <- 1/sd.alpha0^2

alpha1 ~ dt(0, 1/2.5^2, 1)           # parameter for coefficient of substrate

mu.alpha2 ~ dt(0, 1/2.5^2, 1)        # Mean hyperparameter for coefficient of altitude.st
sd.alpha2 ~ dt(0, 1/25^2, 1) T(0,) # SD hyperparameter
tau.alpha2 <- 1/sd.alpha2^2

mu.alpha3 ~ dt(0, 1/2.5^2, 1)        # Mean hyperparameter for coefficient of altitude.st .2
sd.alpha3 ~ dt(0, 1/25^2, 1) T(0,) # SD hyperparameter
tau.alpha3 <- 1/sd.alpha3^2

mu.alpha4 ~ dt(0, 1/2.5^2, 1)        # Mean hyperparameter for coefficient of precipiation.st
sd.alpha4 ~ dt(0, 1/25^2, 1) T(0,) # SD hyperparameter
tau.alpha4 <- 1/sd.alpha4^2

mu.alpha5 ~ dt(0, 1/2.5^2, 1)        # Mean hyperparameter for coefficient of precipitation.st.2
sd.alpha5 ~ dt(0, 1/25^2, 1) T(0,) # SD hyperparameter
tau.alpha5 <- 1/sd.alpha5^2

mu.alpha6 ~ dt(0, 1/2.5^2, 1)        # Mean hyperparameter for coefficient of solar.st
sd.alpha6 ~ dt(0, 1/25^2, 1) T(0,) # SD hyperparameter
tau.alpha6 <- 1/sd.alpha6^2

mu.alpha7 ~ dt(0, 1/2.5^2, 1)        # Mean hyperparameter for coefficient of solar.st.2
sd.alpha7 ~ dt(0, 1/25^2, 1) T(0,) # SD hyperparameter
tau.alpha7 <- 1/sd.alpha7^2



## Priors for detection process

for(o in 1:nobservers){
  mu.beta1.prior[o] ~ dunif(0,1)         # Mean hyperparameters for observers
  mu.beta1[o] <- logit(mu.beta1.prior[o])
  sd.beta1[o] ~ dt(0, 1/25^2, 1) T(0,) # SD hyperparameters
  tau.beta1[o] <- 1/sd.beta1[o]^2
}

beta2 ~ dt(0, 1/2.5^2, 1)     # prior for coefficient of conspicuousness
beta3 ~ dt(0, 1/2.5^2, 1)     # prior for coefficient of identifiability
beta4 ~ dt(0, 1/2.5^2, 1)     # prior for coefficient of experience

}
")


# Initial values
zst <- matrix(NA, nsites, nspecies)
for(i in 1:nspecies){
  zst[,i] <- apply(data$y[,,i], 1, max, na.rm=TRUE)
}
inits <- function(){list(z=zst, alpha0=rnorm(nspecies),
                         mu.alpha0.prior=runif(1), sd.alpha0=runif(1,0.1,2),
                         alpha1=runif(1,0,2),
                         alpha2=rnorm(nspecies), mu.alpha1=rnorm(1), sd.alpha1=runif(1,0.1,2),
                         alpha3=rnorm(nspecies), mu.alpha2=rnorm(1), sd.alpha2=runif(1,0.1,2),
                         alpha4=rnorm(nspecies), mu.alpha4=rnorm(1), sd.alpha4=runif(1,0.1,2),
                         alpha5=rnorm(nspecies), mu.alpha5=rnorm(1), sd.alpha5=runif(1,0.1,2),
                         alpha6=rnorm(nspecies), mu.alpha6=rnorm(1), sd.alpha6=runif(1,0.1,2),
                         alpha7=rnorm(nspecies), mu.alpha7=rnorm(1), sd.alpha7=runif(1,0.1,2),
                         mu.beta1.prior=runif(nobservers), sd.beta1=runif(nobservers,0.1,2),
                         beta2=runif(1,0,1),
                         beta3=runif(1,0,1),
                         beta4=runif(1,0,1)
                         )}

# Parameters monitored
params <- c("mu.alpha0","sd.alpha0",
            "alpha1",
            "mu.alpha2","sd.alpha2",
            "mu.alpha3","sd.alpha3",
            "mu.alpha4","sd.alpha4",
            "mu.alpha5","sd.alpha5",
            "mu.alpha6","sd.alpha6",
            "mu.alpha7","sd.alpha7",
            "mu.beta1","sd.beta1",
            "beta2",
            "beta3",
            "beta4"
            )

# MCMC settings
na=2 ;  nc=1 ;  ni=4 ;  nb=2 ;  nt=1 # try-out: 20 min
na=6000 ;  nc=5 ;  ni=50000 ;  nb=30000 ;  nt=20 # takes approximately 3 days!
paste(nc*(ni-nb)/nt, "posterior samples")


# Call JAGS
output <- jags(data, inits, params, "1_code/jags_model_priors2.txt", n.chains=nc,
               n.thin=nt, n.iter=ni, n.burnin=nb, parallel=ifelse(nc>1,TRUE,FALSE))
save(output, file="2_output/output.jags_model_priors2.RData")
