#'Logistic Growth Model.
#'@param K - Carrying Capacity
#'@param r - Growth Rate
#'@param N0 - Population Start Size
#'@param t - Time points for which to do inference; single value or vector
#'@return A vector of cell count(s) for the specified input time(s) is returned.
#'@export
detLog=function(K,r,N0,t){
  return(K*N0*exp(r*t)/(K+N0*(exp(r*t)-1)))
}

#'Hybrid model expressed as number of cells at time t1, after starting at t0.
#'@param K - Carrying Capacity
#'@param r - Growth Rate
#'@param N0 - Population Start Size
#'@param Nswitch - Switch from stochastic to deterministic Model
#'@param t0 - Initial time point
#'@param t1 - Time points for which to do inference; single value or vector
#'@return A vector of cell count(s) for the specified input time(s) is returned.
#'@export
simDt=function(K=1000,r=1,N0=1,NSwitch=100,t0=0,t1=1){
  if(NSwitch>N0){
    # Unusually, for this model, we know the number of events a priori
    eventNo=NSwitch-N0
    # So we can just generate all required random numbers (quickly) in one go
    unifs=runif(eventNo)
    clist=seq(N0,NSwitch)
    # Time between events
    dts=-log(1-unifs)/(r*clist[seq(2,(eventNo+1))]*(1-clist[seq(2,(eventNo+1))]/K))
    # Absolute times
    ats=c(t0,t0+cumsum(dts))
    tmax=max(ats)
    if(tmax>=t1){
      # Interpolate for estimate of c at t1
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


#'Bayesian Inference for Hybrid Model.
#'@param params - Initial guess for Growth Rate r and Carrying Capacity k in a vector
#'@param pmin - Lower bound of prior for r and K
#'@param pmax - Upper bound for prior for r and K
#'@param tune - Tuning Parameter
#'@param thin - Thinning
#'@param noiseSD - Precision Measures (Standard Deviation)
#'@param switchN - Number of cell after which to switch from the stochastic to the deterministic model
#'@param dat - Growth curve on which to no inference. Vector of cell counts with time points specified as rownames
#'@return MCMC chain for r and K.
#'@export
BayesHybrid<-function(params,pmin,pmax,iters,tune,thin,noiseSD,switchN,dat){
  noiseSD=noiseSD
  stepSim=function(x0=1, t0=0, deltat=1, th = c(100,3))  simDt(th[1],th[2],x0,switchN,t0,t0+deltat)
  dataLik<-function(x,t,y,log=TRUE,...){
    ll=sum(dnorm(y,x,noiseSD,log=TRUE))
    if(log)
      return(ll)
    else
      return(exp(ll))
  }
  simx0=function(N,t0,...){
    # Initial condition for all populations is 1 cell
    return(matrix(1,nrow=N))
  }
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
          w[j]=dataLik(xmat[j,],t=times[i+1],y=data[i,],log=FALSE,...)
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
  mLLik=pfMLLik(20,simx0,0,stepSim,dataLik,dat)
  # MCMC algorithm
  p=length(params)
  ll=-1e99
  parmat=matrix(0,nrow=iters,ncol=p)
  colnames(parmat)=names(params)
  # Main MCMC loop
  for(i in 1:iters){
    if (i%%(iters/10)==0) message(paste(i,date()))
    for(j in 1:thin){
      paramsprop=pmin/2
      while(sum((paramsprop<pmin)|(paramsprop>pmax))>0){
        # Reject particles outside of range
        paramsprop=params*exp(rnorm(p,0,tune))
      }
      llprop=mLLik(paramsprop)
      if (log(runif(1)) < llprop-ll){
        params=paramsprop
        ll=llprop
      }
    }
    parmat[i,]=params
  }
  return(parmat)
}

#'Plotting the Posterior Probabilities for an MCMC chain.
#'@param pmin - Lower bound of prior for r and K
#'@param pmax - Upper bound for prior for r and K
#'@param parmat - MCMC chain returned from \code{BayesHybrid}
#'@return Plots the posterior probability distributions.
#'@export
PlotPosteriorProb<-function(pmin,pmax,parmat){
  op=par(mfrow=c(2,1))
  curve(dunif(x,pmin[2],pmax[2]),from=pmin[2],to=pmax[2],
        xlab="r (1/h)",ylab="Density",main="Growth Rate",cex.lab=1.5,
        ylim=c(0,max(density(parmat[,2])$y)),lty=2)
  points(density(parmat[,2]),main="",lwd=3,type='l')
  curve(dunif(x,pmin[1],pmax[1]),from=pmin[1],to=pmax[1],
        main="Carrying Capacity",xlab="K (cells)",ylab="Density",cex.lab=1.5,
        ylim=c(0,max(density(parmat[,1])$y)),lty=2)
  points(density(parmat[,1]),main="",lwd=3,type='l')
  legend("topright",legend=c("Prior","Posterior"),lwd=c(1,3),col=c("black","black"),lty=c(2,1),cex=0.8)
  par(op)
}

#'Function to simulate growth curves according to the Hybrid Model
#'@param K - Carrying Capacity
#'@param r - Growth Rate
#'@param N0 - Population Start Size
#'@param Nswitch - Switch from stochastic to deterministic model
#'@param detpts - Time interval used for the deterministic model
#'@export
simCellsHybrid=function(K,r,N0,NSwitch,detpts=100){
  # Every event produces one cell and consumes one unit of nutrients
  if(NSwitch>N0){
    # Unusually, for this model, we know the number of events a priori
    eventNo=NSwitch-N0
    # So we can just generate all required random numbers (quickly) in one go
    unifs=runif(eventNo)
    clist=(N0+1):NSwitch
    # Time between events
    dts=-log(1-unifs)/(r*clist*(1-clist/K))
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

#'Function to simulate growth curves stochastically
#'@param K - Carrying Capacity
#'@param r - Growth Rate
#'@param N0 - Population Start Size
#'@export
simCellsStoch=function(K,r,N0){
  # Unusually, for this model, we know the number of events a priori
  eventNo=K-N0
  # So we can just generate all required random numbers (quickly) in one go
  unifs=runif(eventNo)
  # Every event produces one cell and consumes one unit of nutrients
  clist=(N0+1):K
  nlist=K+1-clist
  # Simulate time between events
  dts=-log(unifs)/(r*clist*(1-clist/K))
  return(data.frame(t=c(0,cumsum(dts)),c=c(N0,clist)))
}

#'Plotting the Posterior Predictives overlayed on the Growth Curve.
#'@param data - The growth curve data as supplied to \code{BayesHybrid}
#'@param parmat - MCMC chain returned from \code{BayesHybrid}
#'@param switchN - Number of cell after which to switch from the stochastic to the deterministic model
#'@param guessK - inital guess for K
#'@return Plots the original growth curves, the posterior predictives and the posterior predictives overlayed on the data.
#'@export
PlotPosteriorPredictive<-function(data,parmat,switchN,guessK){
  op=par(mfrow=c(3,1))
  plot(rownames(data),data$c,xlab="Time (h)",ylab="Cell number",cex.lab=1.5,main=paste("Growth Curve"),
       type='p',pch=16,col="red")
  # Posterior Predictive
  plot(NULL,ylim=c(0,guessK+10000),xlim=c(0,100),xlab="Time (h)",ylab="Cell count", main="Posterior Predictives",cex.lab=1.5)
  for (i in 1:dim(parmat)[1]){
    pospred=simCellsHybrid(parmat[i,1],parmat[i,2],1,switchN)
    lines(pospred$t,pospred$c,col=adjustcolor("red",0.5))
  }
  plot(rownames(data),data$c,xlab="Time (h)",ylab="Cell count", main="Posterior Predictives",
       cex.lab=1.5,type='l')
  for (i in 1:dim(parmat)[1]){
    pospred=simCellsHybrid(parmat[i,1],parmat[i,2],1,switchN)
    lines(pospred$t,pospred$c,col=adjustcolor("red",0.05))
  }
  lines(rownames(data),data$c,type='p',lwd=2,pch=4)
  lines(rownames(data),data$c,type='l')
  par(op)
}
