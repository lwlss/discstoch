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
#library(detstocgrowth)
#library(Rcpp)
library(smfsb)
#library(microbenchmark)

#setwd("~/BayesianInference")
#Rcpp::sourceCpp('~/StochasticFunctions.cpp')

####################################### Functions ########################################################

detLog=function(K,r,c0,t,B=1){
  return(B*K*c0*exp(r*t)/(K+c0*(exp(r*t)-1)))
}

# Hybrid model expressed as number of cells at time t1, after starting at t0
simDt=function(K=1000,r=1,B=1,N0=1,NSwitch=100,t0=0,t1=1){
  if(NSwitch>N0){
    if(B==0){
      return(0)
    }else{
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
        return(detLog(K,r,NSwitch,t1-tmax,B))
      }
    }
  }else{
    return(detLog(K,r,N0,t1-t0,B))
  }
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

# # Marginal likelihood
# pfMLLik_cpp=function (n, simx0, t0, stepFun, dataLik, data)
# {
#   times = c(t0, as.numeric(rownames(data)))
#   deltas = diff(times)
#   area=data$c
#   return(function(...) {
#     xmat = simx0(n, t0, ...)
#     w=matrix(nrow=n)
#     ll = 0
#     for (i in 1:length(deltas)) {
#       # Replace apply function with for loop to avoid vectorising vectorised function
#       # if statements, seq function and a:b notation don't play nice with vectorisation...
#       rcpp_op=pf_cpp(simx0,n,0,times,deltas,dataLik,stepSim,i,area)
#       xmat=rcpp_op$xmat
#       w=rcpp_op$w
#       if (max(w) < 1e-20) {
#         warning("Particle filter bombed") #WHAT DOES THIS MEAN?
#         return(-1e+99)
#       }
#       ll = ll + log(mean(w))
#       rows = sample(1:n, n, replace = TRUE, prob = w)
#       # Typecast to matrix here, as otherwise 1D matrix gets converted to list...
#       xmat = as.matrix(xmat[rows, ])
#     }
#     ll
#   })
# }

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
    #print(i)
    if (i%%(iters/10)==0) message(paste(i,date(),"Bb"))
    for (j in 1:thin) {
      thprob=pmin-1
      while(sum((thprob[1]<=pmin[1])|(thprob[1]>=pmax[1]))>0){
        # Reject particles outside of range
        thprob[1]=th[1]*exp(rnorm(1,0,tune))
      }
      thprob[2]=rbinom(1,1,0.8)
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
    area=fread("~/BayesianInference/Lawless_area_shortTC.txt",header=FALSE)
    times=fread("~/BayesianInference/Lawless_time_shortTC.txt",header=FALSE)
    data=fread("~/BayesianInference/Lawless_data_shortTC.txt",header=FALSE) #3rd column (Identifier) => strain_parentcolony
    names(data)=c("genotype","clonalcolony","identifier","blobnumber")
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=fread("~/BayesianInference/Levy_area_filtered1.txt",header=FALSE)
    times=fread("~/BayesianInference/Levy_times_filtered1.txt",header=FALSE)
    data=fread("~/BayesianInference/Levy_data_filtered1.txt",header=TRUE) #3rd column (Identifier) => replicate
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    area=fread("~/BayesianInference/Ziv_area_filtered1.txt",header=FALSE)
    times=fread("~/BayesianInference/Ziv_times_filtered1.txt",header=FALSE)
    data=fread("~/BayesianInference/Ziv_data_filtered1.txt",header=TRUE) #3rd column (Identifier) => colony
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else {print("Not a valid dataset")}
}

# Choosing a data set
#datsetname="Lawless"
x=dataset(datsetname)
area=x$area
times=x$times
data=x$data
data=as.matrix(data)
area=as.matrix(area)

# Getting the data into the right format
area_cell=median(area[,1]) #Lawless 92.5
#area_cell=16.67 #Levy & Ziv
calibrated_area=t(apply(area,1, function(x) x/area_cell))
#gc=1367
modelled_data=data.frame(c=calibrated_area[gc,],t=t(times[gc,]))
rownames(modelled_data)=times[gc,]
modelled_data$t=NULL
plot(as.numeric(rownames(modelled_data)),modelled_data$c,
     main="Simulated Growth Curve",ylab="Cell Count",xlab="Time (h)",cex.lab=1.4)

#print(detstocgrowth::LM_growthrate(modelled_data$c,rownames(modelled_data))$rate)

# # Plotting all growth curves to choose which ones to model
# pdf(height = 16, width = 16, file = paste("BayesianInferenceGCs_",datsetname,".pdf",sep=""))
# par(mfrow=c(4,4))
# for (i in 1:dim(calibrated_area)[1]){
#   plot(as.numeric(times[i,]),as.numeric(calibrated_area[i,]),
#        #main=paste(i,data[i,1],data[i,2],data[i,3],data[i,4]),ylab="Cell Count",xlab="Time (h)",cex.lab=1.4,type='l',lty=2)
#        main=paste(i,data[i,1],data[i,2],data[i,3]),ylab="Cell Count",xlab="Time (h)",cex.lab=1.4,type='l',lty=2)
#   points(as.numeric(times[i,]),as.numeric(calibrated_area[i,]))
# }
# dev.off()

# Using the hybrid simulation but setting the switch much higher than any of the time coure data
#switchN=65
switchN=round(max(calibrated_area))+1
noiseSD=10

# Step Function for pfMLLik
stepSim=function(x0=1, t0=0, deltat=1, B=1, th=c(3,1))  simDt(K=15000,th[1],th[2],x0,switchN,t0,t0+deltat) #fix K and change dimension of th
#stepSim_cpp=function(x0=1, t0=0, deltat=1, th = c(100,3))  simDt_cpp(th[1],th[2],x0,switchN,t0,t0+deltat)
#is this what is making stepSim_cpp really slow? try putting this is cpp as well....

# Number of particles
n=10

mLLik=pfMLLik(n,simx0,0,stepSim,dataLik,modelled_data)
#mLLik_cpp1=pfMLLik(n,simx0,0,stepSim_cpp,dataLik,modelled_data)
#mLLik_cpp2=pfMLLik_cpp(n,simx0,0,stepSim,dataLik,modelled_data)

# MCMC algorithm
print(date())
# iters=100
# tune=0.01
thin=iters/10
th=c(r=0.2,dfrac=0.8)
p=length(th)
ll=-1e99
#Priors
pmin=c(r=0.1,dfrac=0) 
pmax=c(r=0.6,drafc=1)
# Main pMCMC loop
#thmat=mcmc_cpp(p,tune,iters,thin,mLLik,th,pmin,pmax)
thmat=mcmc(p,tune,iters,thin,mLLik,th,pmin,pmax)
message("Done!")
print(date())
# Compute and plot some basic summaries
print(mcmcSummary(thmat))

pdf(height = 8, width = 9,file = paste(datsetname,"_Bb0106_Stoch_MCMC_Summary_GC_",gc,"_Iters",iters,"_tune",tune,".pdf",sep=""))
#svg(paste("MCMC_Summary_GC",gc,".svg",sep=""),width=7, height=21,pointsize=24)
mcmcSummary(thmat,show=FALSE,plot=TRUE)
op=layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
op=par(mfrow=c(2,1))
plot(as.numeric(times[gc,]),as.numeric(calibrated_area[gc,]),
     main=paste(gc,data[gc,1],data[gc,2],data[gc,3],data[gc,4]),ylab="Cell Count",xlab="Time (h)",cex.lab=1.4,type='l',lty=2)
points(as.numeric(times[gc,]),as.numeric(calibrated_area[gc,]))
curve(dunif(x,pmin[1],pmax[1]),from=pmin[1],to=pmax[1],
      xlab="r (1/h)",ylab="Density",main="Growth Rate",cex.lab=1.5,
      ylim=c(0,max(density(thmat[,1])$y)),lty=2)
points(density(thmat[,1]),main="",lwd=3,type='l')
legend("topright",legend=c("Prior","Posterior"),lwd=c(1,3),col=c("black","black"),lty=c(2,1))
# curve(dunif(x,pmin[2],pmax[2]),from=pmin[1],to=pmax[1],
#       main="Carrying Capacity",xlab="K (cells)",ylab="Density",cex.lab=1.5,ylim=c(0,max(density(thmat[,1])$y)),lty=2)
# points(density(thmat[,1]),main="",lwd=3,type='l')
# legend("topright",legend=c("Prior","Posterior"),lwd=c(1,3),col=c("black","black"),lty=c(2,1))
par(op)
dev.off()

# # Microbenchmarks
# op1=microbenchmark(pfMLLik(10,simx0,0,stepSim,dataLik,modelled_data),
#                      pfMLLik(10,simx0,0,stepSim_cpp,dataLik,modelled_data),
#                      pfMLLik_cpp(10,simx0,0,stepSim,dataLik,modelled_data),
#                      pfMLLik_cpp(10,simx0,0,stepSim_cpp,dataLik,modelled_data))
# op2=microbenchmark(mcmc(p,tune,iters,thin,mLLik,th,pmin,pmax),
#                      mcmc_cpp(p,tune,iters,thin,mLLik,th,pmin,pmax))
#
# print(op1)
# boxplot(op1,names=c("stepSim","stepSim_cpp","pfMLLik", "pfMLLik_cpp"))
# print(op2)
# boxplot(op2,names=c("mcmc","mcmc_cpp"))
