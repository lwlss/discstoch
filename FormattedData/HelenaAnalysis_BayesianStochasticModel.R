# Bayesian Stochastic Model 
setwd("~/BayesianInference")
library(data.table)
library(detstocgrowth)
library(Rcpp)

####################################### Functions ########################################################

detLog=function(K,r,c0,t){
  return(K*c0*exp(r*t)/(K+c0*(exp(r*t)-1)))
}

Rcpp::cppFunction(
  'NumericVector simDt_cpp(int K=1000, double r=1.0, int N0=0, int NSwitch=100, double t0=0, double t1=1) {
  Environment myEnv = Environment::global_env();
  Function detLog = myEnv["detLog"];
    if (NSwitch>N0){
      int eventN0=NSwitch-N0;
      NumericVector unifs=runif(eventN0);
      IntegerVector nn = seq(N0,NSwitch);
      NumericVector clist = as<NumericVector>(nn);
      NumericVector dts=-log(unifs)/(r*clist[seq(2,(eventN0+1))]*(1-clist[seq(2,(eventN0+1))]/K));
      NumericVector ats=cumsum(dts);
      ats.push_front(t0);
      double tmax=max(ats);
      if (tmax>=t1){
        Environment stats("package:stats");
        Function approxfun = stats["approxfun"];
        Rf_PrintValue(ats);
        Rf_PrintValue(clist);
        NumericVector af=approxfun(ats,clist,"constant");
        return af(t1);
      }else{
        return detLog(K,r,NSwitch,t1-tmax);
      }
    }else{
      return detLog(K,r,N0,t1-t0);
    }
  }'
)

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
    if(tmax>=t1){
      # Interpolate for estimate of c at t1
      print(ats)
      print(length(ats))
      print(length(clist))
      af=approxfun(ats,clist,method="constant")
      return(af(t1))	
    }else{
      # Deterministic simulation from tmax to t1
      return(detLog(K,r,NSwitch,t1-tmax))
    }
  }else{
    return(detLog(K,r,N0,t1-t0))
  }
}

# Writing the above function in Rcpp


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
    area=fread("Levy_area_filtered.txt",header=FALSE)
    times=fread("Levy_times_filtered.txt",header=FALSE)
    data=fread("Levy_data_filtered.txt",header=TRUE) #3rd column (Identifier) => replicate
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
datsetname="Lawless"
x=dataset(datsetname)
area=x$area
times=x$times
data=x$data

# Choosing a strain and extracting the data for it
strain_names=unique(data$genotype)
pickstrain=strain_names[1]#choose strain here!
strain=subset_strain(data,area,times,pickstrain)

# Getting the data into the right format
area_cell=17.25
calibrated_area=t(apply(area,1, function(x) x/area_cell))
modelled_data=data.frame(c=calibrated_area[1,],t=t(times[1,]))
rownames(modelled_data)=times[1,]
modelled_data$t=NULL
plot(modelled_data$c,as.numeric(rownames(modelled_data)),
     main="Simulated Growth Curve",ylab="Cell Count",xlab="Time (h)",cex.lab=1.4)

# Using the hybrid simulation but setting the switch much higher than any of the time coure data 
switchN=1000

# Step Function for pfMLLik
stepSim=function(x0=1, t0=0, deltat=1, th = c(100,3))  simDt_cpp(th[1],th[2],x0,switchN,t0,t0+deltat)

# Log likelihood of the observation
dataLik<-function(x,t,y,log=TRUE,...){
  ll=sum(dnorm(y,x,noiseSD,log=TRUE))
  if(log)
    return(ll)
  else
    return(exp(ll))
}

# Prior distribution for X0
simx0=function(N,t0,...){
  # returns a matrix whose rows are random samples from the initial distribution
  return(matrix(sample(calibrated_area[,t0+1],N),nrow=N))
  # or would it make sense to set these all equal to area_cell?
}

# Marginal likelihood
pfMLLik=function (n, simx0, t0, stepFun, dataLik, data)
{
  print("1 Got here fine")
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
        w[j]=dataLik(xmat[j,],t=times[i+1],y=data[i,],log=FALSE,...)
      }
      if (max(w) < 1e-20) {
        warning("Particle filter bombed") #WHAT DOES THIS MEAN?
        return(-1e+99)
      }
      ll = ll + log(mean(w))
      rows = sample(1:n, n, replace = TRUE, prob = w)
      # Typecast to matrix here, as otherwise 1D matrix gets converted to list...
      xmat = as.matrix(xmat[rows, ])
    }
    ll
    print("2 Got here fine")
  })
}

mLLik=pfMLLik(5,simx0,0,stepSim,dataLik,modelled_data)

#Check a best guess in the particle filter
print(mLLik(th=c(1000,0.2))) #K,r

# MCMC algorithm
date()
iters=1000
tune=0.01
thin=10
# Flat priors - change these to uniform priors such that K~(60,1000), r~(0,3)
th=c(K = 1000, r = 0.2)
p=length(th)
ll=-1e99
thmat=matrix(0,nrow=iters,ncol=p)
colnames(thmat)=names(th)
# Main pMCMC loop
for (i in 1:iters) {
  message(paste(i,""),appendLF=FALSE)
  for (j in 1:thin) {
    thprop=th*exp(rnorm(p,0,tune))
    llprop=mLLik(thprop)
    if (log(runif(1)) < llprop - ll) {
      th=thprop
      ll=llprop
    }
  }
  thmat[i,]=th
}
message("Done!")
date()
# Compute and plot some basic summaries
mcmcSummary(thmat)
print(apply(thmat,2,mean))