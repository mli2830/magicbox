library(R2jags)


simdat <- function(seed=1099,beta=0.5,gamma=0.5,s0=49,i0=1,t0=1,
                       tot.t=20,dt=1,N=50,Pc=0.5) {
  if (!is.null(seed)) set.seed(seed)
  tvec <- seq(1,tot.t,by=dt)
  n <- length(tvec)
  S <- numeric(n)
  I <- numeric(n)
  Pi <- numeric(n)
  inf <- numeric(n)
  rec <- numeric(n)
  Iobs <- numeric(n)
  S[1] <- s0
  I[1] <- i0
  Pr <- 1 - exp(-gamma*dt)
  Pi[1] <- 1 - exp(-beta*I[1]*dt/N)
  inf[1] <- rbinom(1,S[1],Pi[1])
  rec[1] <- rbinom(1,I[1],Pr)
  Iobs[1] <- rbinom(1,I[1],Pc)
  for (i in 2:n) {
    S[i] <- S[i-1] - inf[i-1]
    I[i] <- I[i-1] - rec[i-1] + inf[i-1]
    Iobs[i] <- rbinom(1,I[i],Pc)
    Pi[i] <- 1 - exp(-beta*I[i]*dt/N)
    inf[i] <- rbinom(1,S[i],Pi)
    rec[i] <- rbinom(1,I[i],Pr)
  }
  cbind(tvec,S,I,Iobs)
}

sirdat <- simdat(seed=13011,tot.t=20,beta=0.6,gamma=0.3,Pc=0.5)

plot(Iobs~tvec,data=sirdat)


parameters <- c("beta","gamma","Pc")

jagsdata <- list(o=sirdat[,4], dt=1,M=nrow(sirdat),s0=48, i0=2, N=50)

jags.sir.model <- function() {
    ## initial values
    S[1] <- s0
    I[1] <- i0
    Pr <- 1 - exp(-gamma*dt)
    Pi[1] <- 1 - exp(-beta*I[1]*dt/N)
    o[1] ~ dbin(Pc,I[1])
    inf[1] ~ dbin(Pi[1],S[1])
    rec[1] ~ dbin(Pr,I[1])
    ## step through observations
    for (i in 2:M) {
      S[i] <- S[i-1] - inf[i-1]
      I[i] <- I[i-1] - rec[i-1] + inf[i-1]
      if (I[i] < o[i]) {
        I[i] <- o[i]
      }
      Pi[i] <- 1 - exp(-beta*I[i]*dt/N)
      inf[i] ~ dbin(Pi[i],S[i])
      rec[i] ~ dbin(Pr,I[i])
      o[i] ~ dbin(Pc,I[i])
    }
    ## aux variables
    
    ## priors
    beta ~ dunif(0,1)
    gamma ~ dunif(0,1)
    Pc ~ dunif(0,1)
}

inits <- list(list(beta=0.6,gamma=0.3,Pc=0.5))

jagsdata$o

j1 <- jags(data=jagsdata,inits, param=parameters, model.file = jags.sir.model,   
           n.chains=length(inits), n.iter=15000, n.burnin=5000, n.thin=37)



j1



