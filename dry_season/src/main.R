library(parallel)
system.time({
cat(detectCores())      #this line tells you how many cores your computer has.
cluster=makeCluster(spec  = 3, type="PSOCK")
clout=clusterEvalQ(cl=cluster,expr=source("WetMCMCEachChain.R"))      #change based on season analyzed
stopCluster(cl=cluster)
})

#First, we summarize for each cluster. Then we combine data in each clusters and do the summary.
#summary output from cluster 1
n.iter=clout[[1]]$value$n.iter
n.burnin=clout[[1]]$value$n.burnin
sigma=clout[[1]]$value$sigma
lamb0=clout[[1]]$value$lamb0
psi=clout[[1]]$value$psi
N=clout[[1]]$value$N
##trace plot include iterations before burnin
plot(sigma,type="l")
plot(lamb0,type="l")
plot(psi,type="l")
plot(N,type="l")
##mean and quantile
outmatrix1<-cbind(sigma,lamb0,psi,N)
library(coda)
outmcmc1<-mcmc(outmatrix1,start=n.burnin)
summary(outmcmc1)
##plot after burnin
plot(outmcmc1)

#summary output from cluster 2
n.iter=clout[[2]]$value$n.iter
n.burnin=clout[[2]]$value$n.burnin
sigma=clout[[2]]$value$sigma
lamb0=clout[[2]]$value$lamb0
psi=clout[[2]]$value$psi
N=clout[[2]]$value$N
##trace plot include iterations before burnin
plot(sigma,type="l")
plot(lamb0,type="l")
plot(psi,type="l")
plot(N,type="l")
##mean and quantile
outmatrix2<-cbind(sigma,lamb0,psi,N)
library(coda)
outmcmc2<-mcmc(outmatrix2,start=n.burnin)
summary(outmcmc2)
##plot after burnin
plot(outmcmc2)

#summary output from cluster 3
n.iter=clout[[3]]$value$n.iter
n.burnin=clout[[3]]$value$n.burnin
sigma=clout[[3]]$value$sigma
lamb0=clout[[3]]$value$lamb0
psi=clout[[3]]$value$psi
N=clout[[3]]$value$N
##trace plot include iterations before burnin
plot(sigma,type="l")
plot(lamb0,type="l")
plot(psi,type="l")
plot(N,type="l")
##mean and quantile
outmatrix3<-cbind(sigma,lamb0,psi,N)
library(coda)
outmcmc3<-mcmc(outmatrix3,start=n.burnin)
summary(outmcmc3)
##plot after burnin
plot(outmcmc3)

#combine data set and do summary
#when you use n cores to run in parallel, you should combine n data set from n clusters.
sigma=c(clout[[1]]$value$sigma[n.burnin:n.iter],clout[[2]]$value$sigma[n.burnin:n.iter],clout[[3]]$value$sigma[n.burnin:n.iter]) 
lamb0=c(clout[[1]]$value$lamb0[n.burnin:n.iter],clout[[2]]$value$lamb0[n.burnin:n.iter],clout[[3]]$value$lamb0[n.burnin:n.iter])
psi=c(clout[[1]]$value$psi[n.burnin:n.iter],clout[[2]]$value$psi[n.burnin:n.iter],clout[[3]]$value$psi[n.burnin:n.iter])
N=c(clout[[1]]$value$N[n.burnin:n.iter],clout[[2]]$value$N[n.burnin:n.iter],clout[[3]]$value$N[n.burnin:n.iter])
##mean and quantile
outmatrixAll<-cbind(sigma,lamb0,psi,N)
library(coda)
outmcmcAll<-mcmc(outmatrixAll,start=n.burnin)
summary(outmcmcAll)
##plot after burnin
plot(outmcmcAll)
