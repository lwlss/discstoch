# Frequentist Deterministic Models

# Simulating exponential growth

#'Applying the exponential model.
#'
#'@param r - Growth Rate
#'@param t - Time
#'@param N0 - Stating Population (default N0=1)
#'@return A list of population estimates for each time point.
#'@export
simExponential<-function(r,t,N0=1){
  N=N0*exp(r*t)
  return(N)
}

#' Log-linear regression of the exponential model.
#'
#'@param ai - Area vector
#'@param ti - Time vector
#'@return Estimated growth rate (\code{rate}) and
#'  starting population (\code{intercept}) as well as the model fit (\code{fit})
#'  which can be accessed using the \code{$} operator.
#'@export
LM_growthrate<-function(ai,ti){
  ai=as.numeric(ai)[which(!is.na(ai))]
  ti=as.numeric(ti)[which(!is.na(ai))]
  ai=as.numeric(ai)[which(ai>0)]
  ti=as.numeric(ti)[which(ai>0)]
  fit<-lm(log(ai)~ti) #this is applying the exponential model
  rate=fit$coefficient[[2]]
  intercept=fit$coefficient[[1]]
  return(list("rates"=rate,"fit"=fit,"int"=intercept))
}

#' Growth rates estimates inferred using the exponential model on the original scale.
#'
#'@param ai - Area vector where each row consists of one growth curve.
#'@param ti - Time vector where each row corresponds to one growth curve.
#'@return Estimated growth rate (\code{rate}) and
#'  starting population (\code{intercept}).
#'@export
EXP_growthrate<-function(ai,ti){
  ai=as.numeric(ai)[which(!is.na(ai))]
  ti=as.numeric(ti)[which(!is.na(ai))]
  ai=as.numeric(ai)[which(ai>0)]
  ti=as.numeric(ti)[which(ai>0)]
  mod1=nls(y~A*exp(r*x),data=data.frame(x=as.numeric(ti),y=as.numeric(ai)),start=list(A=min(ai),r=0.2))
  k=summary(mod1)$parameters[2,1]
  int=summary(mod1)$parameters[1,1]
  return(list("rate"=k,"int"=int))
}

#' Piece-wise log linear regression of the exponential model.
#'
#'@param ai - Area matrix where each row consists of one growth curve.
#'@param ti - Vector of chosen time intervals for simulating.
#'@return Estimated slope before (\code{slope1}) and after (\code{slope2}) the
#'  breakpoint (\code{bp}) along with the estimated starting population
#'  (\code{intercept}) and fit (\code{fit}) which can be accessed using the \code{$} operator.
#'@export
piecewise_LM<-function(ai,t){
  dat=data.frame(x=t,y=ai)
  lin_mod=lm(y~x,data=dat)
  fit<-segmented::segmented(lin_mod,seg.Z=~x,psi=10,data=dat,
                            model=TRUE)
  #control=seg.control(display=FALSE)
  pb=fit$psi[2]
  #converting intercept values from log(area) to area
  intercept=exp(fit$model$y[which(fit$model$x==round(pb))])
  slope1=segmented::slope(fit)$x[1]
  slope2=segmented::slope(fit)$x[2]
  return(list("slope1"=slope1,"slope2"=slope2,"pb"=pb,"fit"=fit, "intercept"=intercept))
}

#' Overlaying Growth Rate Density Distributions.
#'
#'@param rates - Vector with all growth rates for that data sets
#'@param data - Data for that dataset; dataframe where the first column specifies the genotype
#'@param maxr - Maximum growth rate value displayed on the y-axis
#'@param datsetname - Data set for which the growth rates are specified
#'@return Plot the growth rate density distributions for a dataset.
#'@export
DensityOverlay<-function(rates,data,maxr,datsetname){
  colours=rev(rainbow(length(unique(data$genotype))))
  plot(NULL,xlim=c(0,maxr),ylim=c(0,11),
       xlab="Growth Rate (1/h)", ylab="Density", cex.lab=1.4,
       main=" ")
  for (i in 1:length(unique(data$genotype))){
    strain_rates=rates[which(data$genotype==unique(data$genotype)[i])]
    lines(density(strain_rates,from=0),col=colours[i],lwd=4)
    polygon(density(strain_rates)$x,density(strain_rates)$y,col=adjustcolor(colours[i],0.2))
  }
  if(datsetname=="Levy"){
    leg=c()
    for (i in LevyP){
      leg=c(leg,paste(i,"; N=",length(which(data$genotype==i)),sep=""))
    }
  }else{
    leg=c()
    for (i in unique(data$genotype)){
      leg=c(leg,paste(i,"; N=",length(which(data$genotype==i)),sep=""))
    }
  }
  legend("topright",legend=leg,col=colours,
         lty=rep(1,(length(unique(data$genotype)))),
         lwd=rep(3,(length(unique(data$genotype)))))
}
