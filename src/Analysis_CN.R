############################
## Read and subset the data
############################
  setwd("C:\\Users\\Nadeau\\Google Drive\\DissertationResearch\\ExperimentalEvolution\\CommunityMonopolizationPaper")
  d = read.csv("hfx_monopolization_data_12June2017.csv")
  
  ##keep only the 48 temperature because this is where we expect monopolization
    d.sub = d[d$temp == 48 & d$tmnt %in% c("C", "T"), ]
    d.sub$tmnt = factor(d.sub$tmnt, levels = c("C", "T"))
    boxplot(rel.a ~ tmnt, data = d.sub)

###############################################
## Prepare the data for JAGS and run the model
###############################################
  jags.data = list(y = d.sub$rel.a,
                   T48 = as.numeric(d.sub$tmnt)-1,
                   clone = d.sub$clone,
                   N = nrow(d.sub),
                   C = length(levels(d.sub$clone)))
  
  parms = c("B0", "B1", "phi", "tau.c", "y.pred")
  
  library(R2jags)
  jags.out = jags.parallel(data = jags.data,
                           parameters.to.save = parms,
                           model.file = "JAGS_SimplestModel.R",
                           n.chains = 3,
                           n.iter = 10000,
                           n.burnin = 5000,
                           n.thin = 5) 
  
  traceplot(jags.out)
  jags.out

###############################
## Posterior predictive checks
###############################
  hist(jags.out$BUGSoutput$sims.list$y.pred, col = rgb(0,0,1,0.5), freq = F)
  hist(jags.data$y, col = rgb(1,0,0,0.5), add = T, freq = F)
  
  pred.sd = apply(jags.out$BUGSoutput$sims.list$y.pred, 1, sd)
  hist(pred.sd, col = "lightblue")
  abline(v = sd(jags.data$y), col = "red", lwd = 3)
  
  pred.mean = apply(jags.out$BUGSoutput$sims.list$y.pred, 1, mean)
  hist(pred.mean, col = "lightblue")
  abline(v = mean(jags.data$y), col = "red", lwd = 3)

####################
## Plot the results
####################
  T42.q = quantile(plogis(jags.out$BUGSoutput$sims.list$B0),
                    probs = c(0.025, 0.5, 0.975))

  T48.q = quantile(plogis(jags.out$BUGSoutput$sims.list$B0 +
                          jags.out$BUGSoutput$sims.list$B1),
                   probs = c(0.025, 0.5, 0.975))

  q = rbind(T42.q, T48.q)
  rownames(q) = NULL
  
  bp = barplot(q[,2],
               ylim = c(0, 0.85),
               axes = F,
               ylab = "Proportion H. mediterranei",
               xlab = "Adapted temperautre of H. volcanii")
  axis(side = 2)
  axis(side = 1, at = bp, labels = c("42", "48"))
  arrows(x0 = bp,
         y0 = q[,1],
         y1 = q[,3],
         angle = 90,
         code = 3)  
  