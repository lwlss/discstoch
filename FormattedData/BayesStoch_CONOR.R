
gc=1
iters=1000
tune=0.01

library(data.table)
library(smfsb)

####################################### Functions ########################################################

simExpDt<-function(r,N0=1,t0=0,t1=1){
  return(N0*exp(r*(t1-t0)))
}

# Log likelihood of the observation
dataLik<-function(x,t,y,log=FALSE,...){
  ll=sum(dnorm(y,x,noiseSD,log=FALSE))
  if(log)
    return(ll)
  else
    return(exp(ll))
}

# Marginal likelihood
pfMLLik=function (n, t0, stepFun, dataLik, data)
{
  times = c(t0, as.numeric(rownames(data)))
  deltas = diff(times)
  return(function(...) {
    xmat = 1
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
    #print(i)
    #message(paste(i,""),appendLF=FALSE)
    if (i%%(iters/10)==0) message(paste(i,date(),"E"))
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

# Choosing a data set
area=as.matrix(fread("~/BayesianInference/Ziv_area_filtered1.txt",header=FALSE))
times=fread("~/BayesianInference/Ziv_times_filtered1.txt",header=FALSE)
data=as.matrix(fread("~/BayesianInference/Ziv_data_filtered1.txt",header=TRUE)) #3rd column (Identifier) => colony

# Getting the data into the right format
area_cell=16.67 #Levy & Ziv
calibrated_area=t(apply(area,1, function(x) x/area_cell))
modelled_data=data.frame(c=calibrated_area[gc,],t=t(times[gc,]))
if(sum(is.na(modelled_data))>0){modelled_data=modelled_data[-which(is.na(modelled_data)),]}
rownames(modelled_data)=modelled_data$t
modelled_data$t=NULL
plot(as.numeric(rownames(modelled_data)),modelled_data$c,
     main="Simulated Growth Curve",ylab="Cell Count",xlab="Time (h)",cex.lab=1.4)

# How relevant is this noiseSD value?
noiseSD=10

# Step Function for pfMLLik
stepSim=function(x0=1, t0=0, deltat=1, th = c(3))  simExpDt(th[1],x0,t0,t0+deltat) 

# Number of particles 
n=10

mLLik=pfMLLik(n,0,stepSim,dataLik,modelled_data)

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

pdf(height = 8, width = 9,file = paste(datsetname,"_E01_Stoch_MCMC_Summary_GC_",gc,"_Iters",iters,"_tune",tune,".pdf",sep=""))
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