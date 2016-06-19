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

# Calculating the growth rate (LS for exponential model)

#' Log-linear regression of the exponential model.
#'
#'@param ai - Area matrix where each row consists of one growth curve.
#'@param ti - Time matrix where each row corresponds to one growth curve.
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

# Fitting a piece-wise linear regression using the segmented package

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
