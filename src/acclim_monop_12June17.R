#############################################################################
# R script to analyze Haloferax monopolization pilot experiment conducted on 
# 12 June 2017 by T. Farkas
#############################################################################

# load data community data

md <- read_xlsx(path="~/Dropbox/Projects/Haloferax/Mediterannei_Invasion/Acclim Monop 12June/acclim_monop_12June17_plated_29June.xlsx", sheet = 1)


#md <- md[md$H646e3<100,]
# md <- md[-c(13:20),] # remove failed platings
# md$rel.a <- md$vol_e.6/rowSums(md[,c(5,8)]) # both e-6
# md$rel.a <- md$vol_e.6/(md$vol_e.6 + (.1 * md$med_e.5)) # both with max counts
# md$rel.a.st1 <- md$rel.a / md$innoc_plt1
# md$rel.a.st2 <- md$rel.a / md$innoc_plt2

# # load population data
# 
# md[order(paste(md$clone, md$tmnt)),]
# pd <- read.csv("~/Dropbox/Projects/Haloferax/Mediterannei_Invasion/Acclim Monop 8May/acclim_monop_pops_8May17.csv")
# pd$rel.pale <- pd$cfu_pale_plt2_e7 / pd$cfu_plt2_e7
# means <- aggregate(pd$rel.pale, by=list(pd$tmnt), FUN=mean)
# sds <- aggregate(pd$rel.pale, by=list(pd$tmnt), FUN=function(x) sd(x)/sqrt(length(x)))
# vals <- data.frame(tmnt=means[,1], mean=means[,2], se=sds[,2])
# exs <- c(1,2,4,5)
# vals <- data.frame(vals, xs=exs)
# vals$lower <- vals$mean - vals$se
# vals$upper <- vals$mean + vals$se
# 
# means.mat <- matrix(vals$mean, nrow=2)
# ses.mat <- matrix(vals$se, nrow=2)
# plot results for e-6 

##### H. vol abundance
means <- aggregate(md$H98e7, by=list(md$tt), FUN=mean)
sds <- aggregate(md$H98e7, by=list(md$tt), FUN=function(x) sd(x)/sqrt(length(x)))
vals <- data.frame(tmnt=means[,1], mean=means[,2], se=sds[,2])
exs <- c(1,2,4,5,7,8)

vals <- data.frame(vals, xs=exs)
vals$lower <- vals$mean - vals$se
vals$upper <- vals$mean + vals$se

means.mat <- matrix(vals$mean, nrow=2)
ses.mat <- matrix(vals$se, nrow=2)

#png(file="~/Dropbox/Projects/Haloferax/Mediterannei_Invasion/pilot_monop_results_18APR17.png")
par(mar=c(4, 6, 3, 1))
bp <- barplot(means.mat, beside=TRUE, col=rep(c("lightblue", "pink"), 3),
              #ylab="relative abundance of H. volcanii", 
              ylab="",
              ylim=c(0, 70),
              names.arg=c("C", "F", "T"), cex.lab=1.5, 
              cex.names=1.5, las=1, cex.axis=1.2)
title(ylab="abundance of H. volcanii", line=4.5, cex.lab=1.5)
arrows(x0 = bp, x1=bp, y0=means.mat-ses.mat, y1=means.mat+ses.mat, length=0, 
       code = 3, lwd=2)
legend(x=5, y=70, legend=c("42C", "48C"), fill=c("lightblue", "pink"), 
       bty="n", cex=1, y.intersp=1.2)

# stats

hist(md$H98e7[md$tt=="T42"])

an1 <- glm(H98e7 ~ tmnt * temp, data=md[md$tmnt!="F",], 
           family="poisson")
summary(an1)

an2 <- glm(H98e7 ~ tmnt * temp, data=md[md$tmnt!="F",], 
           family="quasipoisson")
summary(an2)

an2 <- lm(H98e7 ~ tmnt * temp, data=md[md$tmnt!="F",])
summary(an2)
##### relative abundance

md2 <- md[md$H646e3<100,]
md2 <- md
means <- aggregate(md2$rel.a, by=list(md2$tt), FUN=mean)
sds <- aggregate(md2$rel.a, by=list(md2$tt), FUN=function(x) sd(x)/sqrt(length(x)))
vals <- data.frame(tmnt=means[,1], mean=means[,2], se=sds[,2])
exs <- c(1,2,4,5,7,8)

vals <- data.frame(vals, xs=exs)
vals$lower <- vals$mean - vals$se
vals$upper <- vals$mean + vals$se

means.mat <- matrix(vals$mean, nrow=2)
ses.mat <- matrix(vals$se, nrow=2)

# stats

an2 <- lm(rel.a ~ tmnt * temp, data=md2[md2$tmnt!="F",])
summary(an1)

#png(file="~/Dropbox/Projects/Haloferax/Mediterannei_Invasion/pilot_monop_results_18APR17.png")
par(mar=c(4, 7, 3, 1))
bp <- barplot(means.mat, beside=TRUE, col=rep(c("lightblue", "pink"), 3),
              #ylab="relative abundance of H. volcanii", 
              ylab="",
              ylim=c(0, .0001),
              names.arg=c("C", "F", "T"), cex.lab=1.5, 
              cex.names=1.5, las=1, cex.axis=1.2)
title(ylab="relative abundance of H. volcanii", line=5, cex.lab=1.5)
arrows(x0 = bp, x1=bp, y0=means.mat-ses.mat, y1=means.mat+ses.mat, length=0, 
       code = 3, lwd=2)
legend(x=5, y=.0001, legend=c("42C", "48C"), fill=c("lightblue", "pink"), 
       bty="n", cex=1, y.intersp=1.2)

plot(rel.a ~ innoc_plt2, data=md)
summary(lm(rel.a.st2 ~ temp * tmnt, data=md))

# statistical analysis
contrasts(md$temp) <- c(0,1) # focus on 42C
an1 <- lm(rel.a ~ temp * tmnt + innoc, data=md)
summary(an1)

ra <- cbind(md$vol_e.6, (1 * md$med_e.6))
an2 <- glm(ra ~ temp * tmnt + innoc, data=md, family=quasibinomial)
summary(an2)

#

##### stuff related to pilot monop test (18 April) all useless here! ####
md <- md[rowSums(md[,5:8], na.rm=TRUE)>0,]
dils <- c(1, 10, 100, 1000, 10, 100, 1000) # dilution factors 

# get table of abundances for each dilution
abs <- t(apply(md[,5:ncol(md)], MARGIN = 1, function(x) dils*x))
colnames(abs) <- sapply(colnames(abs), function(x) paste("a.", x, sep=""))

md1 <- cbind(md, abs)

lag <- function(r1) {
  
  r1 <- t(r1)  
  lag <- vector()
  for (i in c(1:7)) {
    ell <- r1[i+1]/r1[i]
    ell <- ifelse(is.nan(ell), NA, ell)
    lag[i+1] <- ell
  }
  lag <- lag[c(1:7)]
  lag[5] <- NA
  return(lag)
  
}

dils <- cbind(md[,1:4], t(apply(md[,5:ncol(md)], MARGIN=1, lag)))
apply(dils[,5:11], MARGIN=2, FUN=function(x) aggregate(x, by=list(dils[,4]), FUN=function(x) mean(x, na.rm = TRUE)))

apply(md1[,5:8], MARGIN=1, function(x) max(x, na.rm=TRUE))

md_v <- md[,1:8]
md_m <- md[,c(1:4, 9:11)]

col_ind.v <- function(x) which(x == max(x[5:8], na.rm=TRUE))[1]
col_ind.m <- function(x) which(x == max(x[5:7], na.rm=TRUE))[1]

inds <- apply(md_v, MARGIN=1, col_ind.v)

dills <- as.numeric(substring(colnames(md_v)[inds], 7, 7))-4
dills2 <- 10^(dills)

md_v$col <- as.numeric(inds)
cnt <- apply(md_v, MARGIN=1, function(x) {
  ind <- as.numeric(x[9])
  return(as.numeric(x[ind]))
})

md$vol.a <- cnt*dills2

###
inds <- apply(md_m, MARGIN=1, col_ind.m)

dills <- as.numeric(substring(colnames(md_m)[inds], 7, 7))-4
dills2 <- 10^(dills)

md_m$col <- as.numeric(inds)
cnt <- apply(md_m, MARGIN=1, function(x) {
  ind <- as.numeric(x[8])
  return(as.numeric(x[ind]))
})

md$med.a <- cnt*dills2
md$tot.a <- rowSums(md[12:13])
md$rel.a <- md$vol.a/md$tot.a

means <- aggregate(md$rel.a, by=list(md$tt), FUN=mean)
sds <- aggregate(md$rel.a, by=list(md$tt), FUN=function(x) sd(x)/sqrt(length(x)))
vals <- data.frame(tmnt=means[,1], mean=means[,2], se=sds[,2])
exs <- c(4, 1, 7, 5, 2, 8)
vals <- data.frame(vals, xs=exs)
vals$lower <- vals$mean - vals$se
vals$upper <- vals$mean + vals$se
vals <- vals[order(vals$xs),]

means.mat <- matrix(vals$mean, nrow=2)
ses.mat <- matrix(vals$se, nrow=2)

#png(file="~/Dropbox/Projects/Haloferax/Mediterannei_Invasion/pilot_monop_results_18APR17.png")
par(mar=c(4, 5, 3, 1))
bp <- barplot(means.mat, beside=TRUE, col=rep(c("lightblue", "pink"), 3),
              ylim=c(0, 1), ylab="relative abundance of H. volcanii",
              names.arg=c("founders", "control", "trend"), cex.lab=1.5, 
              cex.names=1.5, las=1, cex.axis=1.2)
arrows(x0 = bp, x1=bp, y0=means.mat-ses.mat, y1=means.mat+ses.mat, length=0, code = 3, lwd=2)
legend(x=bp[4], y=.8, legend=c("42C", "48C"), fill=c("lightblue", "pink"), 
       bty="n", cex=1.5, y.intersp=1.2)
dev.off()

#prediction
pred.mat <- matrix(c(.8, .2, .8, .2, .6, .4), nrow=2)
par(mar=c(4, 5, 3, 1))
bp <- barplot(pred.mat, beside=TRUE, col=rep(c("lightblue", "pink"), 3),
              ylim=c(0, 1), ylab="relative abundance of H. volcanii",
              names.arg=c("founders", "control", "trend"), cex.lab=1.5, 
              cex.names=1.5, las=1, cex.axis=1.2)
legend(x=bp[4], y=1, legend=c("42C", "48C"), fill=c("lightblue", "pink"), 
       bty="n", cex=1.5, y.intersp=1.2)
# prediction

##### analysis #####
md1 <- md[md$tmnt != "found",]
ra <- cbind(md1$vol.a, md1$med.a)
an1 <- lm(rel.a ~ temp * tmnt, data=md1)
summary(an1)
