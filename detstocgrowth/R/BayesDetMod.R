#'Deterministic Bayesian inference on lineage growth curves.
#'
#'@param area - Area Estimates on which to do inference in a matrix
#'@param times - Times associated with Area Observations in a matrix
#'@param model_specification - Model type. 1=logistic; 2=exponential;
#'3=log-linear; 4=exponential with B; 5=log-linear with B
#'@param it - Number of Iterations
#'@param th - Thinning
#'@param ch - Number of Chains
#'@return Returns the mean parameter estimates for all input time courses and the MCMC chain.
#'@export
BayesDet<-function(area,times,model_specification,it,th,ch){
  # Logistic Model
  modelstring1="
  model{
  for (i in 1:n) {
  mu[i]<-(K*x0*exp(r*t[i]))/(K+x0*(exp(r*t[i])-1))
  x[i]~dnorm(mu[i],tau)
  }
  for (j in 1:n) {
  pmu[j]<-(K*x0*exp(r*t[j]))/(K+x0*(exp(r*t[j])-1))
  px[j]~dnorm(pmu[j],tau)
  }
  K~dunif(1000,1000000000)
  r~dunif(0,2)
  x0~dunif(4,300)
  tau~dunif(0,1000)
  }
  "
  # Exponential Model
  modelstring2="
  model{
  for (i in 1:n){
  mu[i]=x0*(exp(r*t[i]))
  x[i]~dnorm(mu[i],tau)
  }
  for (j in 1:n) {
  pmu[j]<-x0*(exp(r*t[j]))
  px[j]~dnorm(pmu[j],tau)
  }
  r~dunif(0,2)
  x0~dunif(4,300)
  tau~dunif(0,1000)
  }
  "
  # Log-Linear Model
  modelstring3="
  model{
  for (i in 1:n){
  mu[i]=log(x0)+(r*t[i])
  x[i]~dlnorm(mu[i],tau)
  }
  for (j in 1:n) {
  pmu[j]=log(x0)+(r*t[j])
  px[j]~dlnorm(pmu[j],tau)
  }
  r~dunif(0,2)
  x0~dunif(4,300)
  tau~dunif(0,1000)
  }
  "
  # Exponential with B
  modelstring4="
  model{
  for (i in 1:n){
  mu[i]=x0*(exp(B*r*t[i]))
  x[i]~dnorm(mu[i],tau)
  }
  for (j in 1:n) {
  pmu[j]<-x0*(B*exp(r*t[j]))
  px[j]~dnorm(pmu[j],tau)
  }
  r~dunif(0.1,0.6)
  x0~dunif(4,300)
  tau~dunif(0,1)
  B~dbern(0.8)
  }
  "
  # Log-linear with B
  modelstring5="
  model{
  for (i in 1:n){
  mu[i]=log(x0)+(B*r*t[i])
  x[i]~dlnorm(mu[i],tau)
  }
  for (j in 1:n) {
  pmu[j]=log(x0)+(B*r*t[j])
  px[j]~dlnorm(pmu[j],tau)
  }
  r~dunif(0.1,0.6)
  x0~dunif(4,300)
  tau~dunif(0,1000)
  B~dbern(0.8)
  }
  "
  if (model_specification==1){
    name="BayesDetLogist"
    print(name)
    model=modelstring1
    print("Logistic Model")
  } else if (model_specification==2){
    name="BayesDetExp"
    print(name)
    model=modelstring2
  } else if (model_specification==3){
    name="BayesDetLogLin"
    print(name)
    model=modelstring3
  } else if (model_specification==4){
    name="BayesDetExpB"
    print(name)
    model=modelstring4
  } else if (model_specification==5){
    name="BayesDetLogLinB"
    print(name)
    model=modelstring5
    }else{
    "Wrong model selection!"
    }
  x0_total=c()
  r_total=c()
  k_total=c()
  B_total=c()
  for (i in 1:dim(area)[1]) {
    vals=which(!is.na(area[1,]))
    N=length(vals)
    prednames=sprintf("pmu[%i]",1:N)
    dat=list('x'=as.numeric(area[i,])[vals], 't'=as.numeric(times[i,])[vals], 'n'=N)
    jags<-rjags::jags.model(textConnection(model),
                            data=dat,
                            n.chains=ch)
    update(jags,it)

    if (model_specification==1){
      samples=rjags::coda.samples(model=jags,variable.names=c('r','x0','tau','K', prednames),n.iter=it,thin=th)
      subset=samples[,c("r","K","x0","tau")]
    } else if (model_specification==2|model_specification==3){
      samples=rjags::coda.samples(model=jags,variable.names=c('r','x0','tau', prednames),n.iter=it,thin=th)
      subset=samples[,c("r","x0","tau")]
    } else if (model_specification==4|model_specification==5){
      samples=rjags::coda.samples(model=jags,variable.names=c('r','x0','tau','B', prednames),n.iter=it,thin=th)
      subset=samples[,c("r","x0","tau",'B')]
    }
    plot(subset)
    op=par(mfrow=c(2,2))
    if (model_specification==3|model_specification==5){
      plot(as.numeric(times[i,])[vals],log(as.numeric(area[i,])[vals]),type='l',xlab="Time",ylab="log(Area)",lty=3,
           main=paste("Growth Curve",i),cex.lab=1.5)
      points(as.numeric(times[i,])[vals],log(as.numeric(area[i,])[vals]),type='p',pch=16)
      preds=samples[[1]][,prednames]
      lines(as.numeric(times[i,])[vals],apply(preds,2,mean),type="l",col="red",ylim=c(0,1000),xlab="Time",ylab="Population Size",
            main="Posterior predictive")
      points(as.numeric(times[i,])[vals],apply(preds,2,quantile,0.05),type="l",lty=3,col="red")
      points(as.numeric(times[i,])[vals],apply(preds,2,quantile,0.95),type="l",lty=3,col="red")
    } else{
      plot(as.numeric(times[i,])[vals],as.numeric(area[i,])[vals],type='l',xlab="Time",ylab="log(Area)",lty=3,
           main=paste("Growth Curve",i),cex.lab=1.5)
      points(as.numeric(times[i,])[vals],as.numeric(area[i,])[vals],type='p',pch=16)
      preds=samples[[1]][,prednames]
      lines(as.numeric(times[i,])[vals],apply(preds,2,mean),type="l",col="red",ylim=c(0,1000),xlab="Time",ylab="Population Size",
            main="Posterior predictive")
      points(as.numeric(times[i,])[vals],apply(preds,2,quantile,0.05),type="l",lty=3,col="red")
      points(as.numeric(times[i,])[vals],apply(preds,2,quantile,0.95),type="l",lty=3,col="red")
    }
    #Posterior Distribution for Growth Rate, r
    m_r=summary(samples)[[1]]['r',1]
    print(paste("Mean of r is", m_r))
    curve(dunif(x,0,2),from=0,to=2,main=paste("Growth Rate, r; Mean:",signif(m_r,3)),ylim=c(0,max(density(samples[[1]][,"r"])$y)))
    points(density(samples[[1]][,"r"]),type="l",col="blue")
    abline(v=m_r,col="black",lty=3)
    #Posterior Distribution for Intercept, x0
    m_x0=summary(samples)[[1]]['x0',1]
    print(paste("Mean of x0 is", m_x0))
    curve(dunif(x,4,300),from=4,to=300,main=paste("Intercept, x0; Mean:",signif(m_x0,3)),ylim=c(0,max(density(samples[[1]][,"x0"])$y)))
    points(density(samples[[1]][,"x0"]),type="l",col="blue")
    abline(v=m_x0,col="black",lty=3)
    #Posterior Distribution for Carrying Capacity, K
    if (model_specification==1){
      m_k=summary(samples)[[1]]['K',1]
      print(paste("Mean of k is", m_k))
      curve(dunif(x,1110,70000),from=1110,to=70000,main=paste("Carrying Capacity, K; Mean:",signif(m_k,3)),ylim=c(0,max(density(samples[[1]][,"K"])$y)))
      points(density(samples[[1]][,"K"]),type="l",col="blue")
      abline(v=m_k,col="black",lty=3)
      k_total=k_total=c(k_total,m_k)
    }
    #Posterior Distribution for Dividing Cells, B
    if (model_specification==4|model_specification==5){
      m_B=summary(samples)[[1]]['B',1]
      print(paste("Mean of B is", m_B))
      curve(dunif(x,0,1),from=0,to=1,main=paste("Dividing Cells, B; Mean:",signif(m_B,3)),ylim=c(0,max(density(samples[[1]][,"B"])$y)))
      points(density(samples[[1]][,"B"]),type="l",col="blue")
      abline(v=m_B,col="black",lty=3)
      B_total=k_total=c(B_total,m_B)
    }
    x0_total=c(x0_total,m_x0)
    r_total=c(r_total,m_r)
  }
  if (model_specification==1){
    Bayes_parameters=data.frame("Intercept"=x0_total,"Rate"=r_total,"CarryingCapacity"=k_total)
  } else if (model_specification==2|model_specification==3){
    Bayes_parameters=data.frame("Intercept"=x0_total,"Rate"=r_total)
  } else if (model_specification==4|model_specification==5){
    Bayes_parameters=data.frame("Intercept"=x0_total,"Rate"=r_total,"DividingCells"=B_total)
  }
  return(list("Parameters"=Bayes_parameters,"Samples"=samples))
}
