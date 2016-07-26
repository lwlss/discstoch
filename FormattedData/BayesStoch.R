# Bayesian Stochastic Model 

# Command line arguments
input <- commandArgs(TRUE)
print(input)

# Print help when no arguments
if(length(input) < 1) {
  input <- c("--help")
}

# Help
if("--help" %in% input) {
  cat("
      The R Script

      Arguments:
      --arg1=someName              - name of the dataset
      --arg2=someValue             - which growth curves
      --arg3=someValue             - number of iterations 
      --arg4=someValue             - tuning parameter
      --help                       - prints this text

      Example:
      ./test.R --arg1='Lawless' --arg2=101
      ")
  
  q(save="no")
}

## Parsing input arguments
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
inputDF <- as.data.frame(do.call("rbind", parseArgs(input)))
inputL <- as.list(as.character(inputDF$V2))
names(inputL) <- inputDF$V1
datsetname=inputL$arg1
gc=as.numeric(inputL$arg2)
iters=as.numeric(inputL$arg3)
tune=as.numeric(inputL$arg4)

library(data.table)
library(smfsb)

####################################### Functions ########################################################

detLog=function(K,r,c0,t){
  return(K*c0*exp(r*t)/(K+c0*(exp(r*t)-1)))
}

# Hybrid model expressed as number of cells at time t1, after starting at t0
simDt=function(K=1000,r=1,N0=1,NSwitch=100,t0=0,t1=1){
  if(NSwitch>N0){
    # Unusually, for this model, we know the number of events a priori
    eventNo=NSwitch-N0
    # So we can just generate all required random numbers (quickly) in one go
    unifs=runif(eventNo)
    clist=seq(N0,NSwitch)
    # Time between events
    dts=-log(unifs)/(r*clist[seq(2,(eventNo+1))]*(1-clist[seq(2,(eventNo+1))]/K))
    # Absolute times
    ats=c(t0,t0+cumsum(dts))
    tmax=max(ats)
    if(tmax>=t1){ #assume max t from stochastic is same as t1 or do linear interpolation 
      # Interpolate for estimate of c at t1
      af=approxfun(ats,clist,method="constant")
      return(af(t1))	
    }else{
      # Deterministic simulation from tmax to t1
      return(detLog(K,r,NSwitch,t1-tmax)) #use analytic solution to test 
    }
  }else{
    return(detLog(K,r,N0,t1-t0))
  }
}

simCellsHybrid=function(K,r,N0,NSwitch,detpts=100){
  # Every event produces one cell and consumes one unit of nutrients
  if(NSwitch>N0){
    # Unusually, for this model, we know the number of events a priori
    eventNo=NSwitch-N0
    # So we can just generate all required random numbers (quickly) in one go
    unifs=runif(eventNo)
    clist=(N0+1):NSwitch
    # Time between events
    dts=-log(unifs)/(r*clist*(1-clist/K))
    # Absolute times
    ats=cumsum(dts)
    tmax=max(ats)
  }else{
    clist=c()
    ats=c()
    tmax=0
  }
    
  # Switch to discrete deterministic logistic function
  clistdet=seq(NSwitch+(K-NSwitch)/detpts,K,(K-NSwitch)/detpts)
  tsdet=log((clistdet*(K - NSwitch))/((K - clistdet)*NSwitch))/r
  return(data.frame(t=c(0,ats,tmax+tsdet),c=c(N0,c(clist,clistdet))))
}

# Log likelihood of the observation
dataLik<-function(x,t,y,log=FALSE,...){
  ll=sum(dnorm(y,x,noiseSD,log=FALSE))
  if(log)
    return(ll)
  else
    return(exp(ll))
}

# Prior distribution for X0
simx0=function(N,t0,...){
  # returns a matrix whose rows are random samples from the initial distribution
  return(matrix(rep(1,N),nrow=N)) #round(sample(calibrated_area[,t0+1],N))
  # or would it make sense to set these all equal to area_cell?
}

# Marginal likelihood
pfMLLik=function (n, simx0, t0, stepFun, dataLik, data)
{
  times = c(t0, as.numeric(rownames(data)))
  deltas = diff(times)
  return(function(...) {
    xmat = simx0(n, t0, ...)
    w=matrix(nrow=n)
    ll = 0
    for (i in 1:length(deltas)) {
      # Replace apply function with for loop to avoid vectorising vectorised function
      # if statements, seq function and a:b notation don't play nice with vectorisation...
      for(j in 1:n) {
        xmat[j,]=stepFun(x0=xmat[j,],t0=times[i],deltat=deltas[i],...)
        w[j]=dataLik(xmat[j,],t=times[i+1],y=data[i,],log=FALSE,...) #likelihood
      }
      if (max(w) < 1e-20) {
        warning("Particle filter bombed") 
        return(-1e+99)
      }
      ll = ll + log(mean(w))
      rows = sample(1:n, n, replace = TRUE, prob = w)
      # Typecast to matrix here, as otherwise 1D matrix gets converted to list...
      xmat = as.matrix(xmat[rows, ])
    }
    ll
  })
}

# Main MCMC loop
mcmc = function(p,tune,iters,thin,mLLik,th,pmin,pmax){
  thmat=matrix(0,nrow=iters,ncol=p)
  colnames(thmat)=names(th)
  for (i in 1:iters) {
    #message(paste(i,""),appendLF=FALSE)
    if (i%%(iters/10)==0) message(paste(i,date()))
    for (j in 1:thin) {
      thprob=pmin-1
      while(sum((thprob<pmin)|(thprob>pmax))>0){
        # Reject particles outside of range
        thprob=th*exp(rnorm(p,0,tune)) 
      }
      llprob=mLLik(thprob)
      if (log(runif(1)) < (llprob - ll)){
        th=thprob
        ll=llprob
      }
    }
    thmat[i,]=th
  }
  return(thmat)
}

###################################### Main ##############################################################

dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=fread("Lawless_area_shortTC.txt",header=FALSE)
    times=fread("Lawless_time_shortTC.txt",header=FALSE)
    data=fread("Lawless_data_shortTC.txt",header=FALSE) #3rd column (Identifier) => strain_parentcolony
    names(data)=c("genotype","clonalcolony","identifier","blobnumber")
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=fread("Levy_area_filtered1.txt",header=FALSE)
    times=fread("Levy_times_filtered1.txt",header=FALSE)
    data=fread("Levy_data_filtered1.txt",header=TRUE) #3rd column (Identifier) => replicate
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    area=fread("Ziv_area_filtered1.txt",header=FALSE)
    times=fread("Ziv_times_filtered1.txt",header=FALSE)
    data=fread("Ziv_data_filtered1.txt",header=TRUE) #3rd column (Identifier) => colony
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else {print("Not a valid dataset")}
}

# Choosing a data set
x=dataset(datsetname)
area=x$area
times=x$times
data=x$data
data=as.matrix(data)
area=as.matrix(area)

# Getting the data into the right format
#area_cell=median(area[,1]) #Lawless 92.5
area_cell=16.67 #Levy & Ziv
calibrated_area=t(apply(area,1, function(x) x/area_cell))
modelled_data=data.frame(c=calibrated_area[gc,],t=t(times[gc,]))
if(sum(is.na(modelled_data))>0){modelled_data=modelled_data[-which(is.na(modelled_data)),]}
rownames(modelled_data)=modelled_data$t
modelled_data$t=NULL
plot(as.numeric(rownames(modelled_data)),modelled_data$c,
     main="Simulated Growth Curve",ylab="Cell Count",xlab="Time (h)",cex.lab=1.4)

# Using the hybrid simulation but setting the switch much higher than any of the time coure data 
#switchN=round(max(calibrated_area))+1
switchN=1
noiseSD=10

# Step Function for pfMLLik
stepSim=function(x0=1, t0=0, deltat=1, th = c(3))  simDt(K=15000,th[1],x0,switchN,t0,t0+deltat) #fix K and change dimension of th 

# Number of particles 
n=10

mLLik=pfMLLik(n,simx0,0,stepSim,dataLik,modelled_data)

# MCMC algorithm
print(date())
thin=iters/10
th=c(r=0.4)
p=length(th)
ll=-1e99
#Priors
pmin=c(r=0) 
pmax=c(r=1)
# Main pMCMC loop
thmat=mcmc(p,tune,iters,thin,mLLik,th,pmin,pmax)
message("Done!")
print(date())
# Compute and plot some basic summaries
print(mcmcSummary(thmat,plot=FALSE))

pdf(height = 8, width = 9,file = paste(datsetname,"_000106_Stoch_MCMC_Summary_GC_",gc,"_Iters",iters,"_tune",tune,".pdf",sep=""))
mcmcSummary(thmat,show=FALSE,plot=TRUE)
op=layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
curve(dunif(x,pmin[1],pmax[1]),from=pmin[1],to=pmax[1],
      xlab="r (1/h)",ylab="Density",main="Growth Rate",cex.lab=1.5,
      ylim=c(0,max(density(thmat[,1])$y)),lty=2)
points(density(thmat[,1]),main="",lwd=3,type='l')
legend("topright",legend=c("Prior","Posterior"),lwd=c(1,3),col=c("black","black"),lty=c(2,1))
op=par(mfrow=c(3,1))
# Posterior Predictive
plot(as.numeric(times[gc,]),as.numeric(calibrated_area[gc,]),
     main=paste(gc,data[gc,1],data[gc,2],data[gc,3]),ylab="Cell Count",xlab="Time (h)",cex.lab=1.5,type='l',lty=2)
points(as.numeric(times[gc,]),as.numeric(calibrated_area[gc,]))
plot(NULL,ylim=c(0,15000),xlim=c(0,200),xlab="Time (h)",ylab="Cell count", main="Posterior Predictive",cex.lab=1.5)
for (i in 1:dim(thmat)[1]){
  pospred=simCellsHybrid(15000,thmat[i,1],1,switchN)
  lines(pospred$t,pospred$c,col=adjustcolor("red",0.5))
}
plot(as.numeric(times[gc,]),as.numeric(calibrated_area[gc,]),
     main="Posterior Predictive Overlay",ylab="Cell Count",xlab="Time (h)",cex.lab=1.5,type='l',lty=1, lwd=3)
for (i in 1:dim(thmat)[1]){
  pospred=simCellsHybrid(15000,thmat[i,1],1,switchN)
  lines(pospred$t,pospred$c,col=adjustcolor("red",0.1))
}
par(op)
dev.off()
