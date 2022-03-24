model{
  for(i in 1:N){
    y[i] ~ dbeta(alpha[i], beta[i])
    alpha[i] <- mu[i] * phi
    beta[i] <- (1 - mu[i]) * phi
    logit(mu[i]) <- B0 + B1*T48[i] + B2[clone[i]]
    
    y.pred[i] ~ dbeta(alpha[i], beta[i])
  }

  for(c in 1:C){
    B2[c] ~ dnorm(0, tau.c)
  }
    
  B0 ~ dnorm(0, 0.001)
  B1 ~ dnorm(0, 0.001)
  phi ~ dgamma(0.01, 0.01)
  tau.c ~ dgamma(0.01, 0.01)
}