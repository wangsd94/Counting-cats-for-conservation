#set parameter
#setwd("E:/wisc/research/Leopard/updated10.26")        #set your working directory
#setwd ("C:/Users/Max Allen/Dropbox/1- Max Allen/Desktop/Top Left- Working Papers/Leopard Density/Analyses/Density") #set your working directory #set your working directory
setwd("C:/Users/maxallen/OneDrive - University of Illinois - Urbana/Manuscripts/Working/Leopard Density/Analyses/Density")
#sink(paste0("E:/wisc/research/Leopard/updated10.26/checkIter", Sys.getpid(), ".txt"))
n.iter=100000                                             #number of iterations
n.burnin=20000                                           #set burnin here, then use in the summary
n.trap=115                                            #number of traps
M=350                                                 #maximum number of individuals
K=3                                                   #number of trapping occassions

#set parameter of prior
sigma_mean = 0.53 #set mean of sigma
sigma_sd = 0.1 #set sd of sigma
lamb0_shape = 2  #set parameter of prior of lamb0
lamb0_rate = 2
psi_shape1 = 4  #set parameter of prior of psi
psi_shape2 = 8
xlim = c(0,22)  #set study area
ylim = c(0,22)
#
s<-array(data=NA,dim=c(M,2,n.iter))
lamb0<-numeric(n.iter)
sigma<-numeric(n.iter)
N<-numeric(n.iter)
psi<-numeric(n.iter)
z<-matrix(data=NA,ncol=n.iter,nrow=M)
bigLambda<-numeric(n.trap)
distsq<-matrix(data=NA,nrow=M,ncol=n.trap)
newBigLambda=numeric(n.trap)
delta.sigma = 0.03
delta.s = 0.1


#read data
n <- read.csv("NLeopDry.csv", header=FALSE)   #encounter history data
n <- as.matrix(n)
X <- read.csv("XLeopDry.csv", header=TRUE)            #trap locations
X <- as.matrix(X)
#write helping function with c++
updateDistance<-function(M,J,S,X){
  D2=matrix(nrow = M, ncol = J)
  S1=matrix(data=S[,1],nrow=M,ncol=J,byrow=FALSE)
  S2=matrix(data=S[,2],nrow=M,ncol=J,byrow=FALSE)
  X1=matrix(data=X[,1],nrow=M,ncol=J,byrow=TRUE)
  X2=matrix(data=X[,2],nrow=M,ncol=J,byrow=TRUE)
  D2=(S1-X1)^2+(S2-X2)^2
  return(D2)
}


#set initial values
sigma[1] = rnorm(1,sigma_mean,sigma_sd)
lamb0[1] = rgamma(1,lamb0_shape,lamb0_rate)
psi[1] = rbeta(1,psi_shape1,psi_shape2)
s[,1,1] = runif(M,xlim[1],xlim[2])
s[,2,1] = runif(M,ylim[1],ylim[2])
z[,1] = rbinom(M,1,psi[1])
distsq<-updateDistance(M,n.trap,as.matrix(s[,,1]),X)
bigLambda = colSums(lamb0[1]*exp(-distsq/(2*sigma[1]^2))*z[,1])
loglike = sum(dpois(n,bigLambda,log = TRUE),na.rm = TRUE)

#initiate acceptance rate
sigma.ups = 0
s.ups = 0
z.ups = 0
#test.lamb.param<-array(dim=c(n.trap,n.iter))


#run mcmc
for(iter in 2:n.iter){
  if(iter%%100==0){cat("this is iteration", iter,"\n")}
  #update sigma
  ###set candidate of sigma in next iteration
  sigma.cand = rnorm(1,sigma[iter-1],delta.sigma)
  while(sigma.cand < 0){  #if candidate of sigma go out of bound, sample again
    sigma.cand = rnorm(1,sigma[iter-1],delta.sigma)
  }
  logprior = dnorm(sigma[iter-1],sigma_mean,sigma_sd,log = TRUE)                  #set prior distribution for sigma
  newBigLambda = colSums(lamb0[iter-1]*exp(-distsq/(2*sigma.cand^2))*z[,iter-1])
  loglike.cand = sum(dpois(n,newBigLambda,log = TRUE),na.rm = TRUE)
  logprior.cand = dnorm(sigma.cand,sigma_mean,sigma_sd,log = TRUE)
  if(runif(1)<exp((loglike.cand+logprior.cand)-(loglike+logprior))){
    sigma[iter]=sigma.cand
    sigma.ups=sigma.ups+1
    bigLambda=newBigLambda
    loglike=loglike.cand
  }else{
    sigma[iter]=sigma[iter-1]
  }
  #update lamb0
  ##when we set prior of lamb0 to be gamma distribution, we have simple posterior(gamma).
  C=bigLambda/lamb0[iter-1]
  #test.lamb.param[,iter-1]=C
  lamb0[iter]=rgamma(1,sum(n,na.rm=T)+lamb0_shape,sum(C)*K+lamb0_rate)
  #update psi
  ##when we set prior of psi to be beta distribution, we have simple posterior(beta)
  psi[iter]=rbeta(1,sum(z[,iter-1])+psi_shape1,M-sum(z[,iter-1])+psi_shape2)
  #update s
  ###set candidate of s in next iteration
  #lambda0 may be updated, so we have to update bigLambda and loglike again
  bigLambda = colSums(lamb0[iter]*exp(-distsq/(2*sigma[iter]^2))*z[,iter-1])
  loglike=sum(dpois(n,bigLambda,log = TRUE),na.rm = TRUE)
  ss=s[,,iter-1]
  for(i in 1:M){
    ##attention here: skip when z[i]==0
    if(z[i,iter-1]==0){
      next
    }
    s.cand=ss
    s.cand[i,]=c(rnorm(1,ss[i,1],delta.s),rnorm(1,ss[i,2],delta.s))
    s.cand.in =s.cand[i,1]>xlim[1]&s.cand[i,1]<xlim[2]&s.cand[i,2]>ylim[1]&s.cand[i,2]<ylim[2]
    while(!s.cand.in){
      s.cand[i,]=c(rnorm(1,ss[i,1],delta.s),rnorm(1,ss[i,2],delta.s))
      s.cand.in =s.cand[i,1]>xlim[1]&s.cand[i,1]<xlim[2]&s.cand[i,2]>ylim[1]&s.cand[i,2]<ylim[2]
    }
    tempdistsq<-updateDistance(M,n.trap,as.matrix(s.cand),X)
    newBigLambda = colSums(lamb0[iter]*exp(-tempdistsq/(2*sigma[iter]^2))*z[,iter-1])
    loglike.cand = sum(dpois(n,newBigLambda,log = TRUE),na.rm = TRUE)
    if(runif(1)<exp(loglike.cand-loglike)){
      ss=s.cand
      s.ups=s.ups+1
      bigLambda=newBigLambda
      loglike = loglike.cand
      distsq=tempdistsq
    }
  }
  s[,,iter]=ss
  
  
  #update z
  zz=z[,iter-1]
  for(i in 1:M){
    #canidate of z (update zi in each iteration)
    z.cand = zz
    z.cand[i] = ifelse(zz[i]== 0, 1, 0)
    logprior = dbinom(zz[i],1,psi[iter],log = TRUE)
    newBigLambda = colSums(lamb0[iter]*exp(-distsq/(2*sigma[iter]^2))*z.cand)
    loglike.cand = sum(dpois(n,newBigLambda,log = TRUE),na.rm = TRUE)
    logprior.cand = dbinom(z.cand[i],1, psi[iter],log = TRUE)
    if(runif(1)<exp((loglike.cand+logprior.cand)-(loglike+logprior))){
      zz = z.cand
      z.ups=z.ups+1
      bigLambda=newBigLambda
      loglike = loglike.cand
    }
  }
  z[,iter] = zz
  N[iter] = sum(zz)
}

out=list(sigma=sigma,lamb0=lamb0,psi=psi,s=s,z=z,N=N,n.iter=n.iter,n.burnin=n.burnin,sigma.ups = sigma.ups,s.ups = s.ups,z.ups =z.ups)
