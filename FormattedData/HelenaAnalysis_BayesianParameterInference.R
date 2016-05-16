# Bayesian Parameter Inference
# Set working directory to source file location 
# http://www4.stat.ncsu.edu/~reich/st590/code/regJAGS
# http://www.r-bloggers.com/parse-arguments-of-an-r-script/

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
      --arg2=someValue:someValue   - which growth curves
      --help                       - prints this text

      Example:
      ./test.R --arg1='Lawless' --arg2=1:100
      ")
  
  q(save="no")
}
## Parsing input arguments
parseArgs <- function(x) strsplit(sub("^--", "", x), "=") 
inputDF <- as.data.frame(do.call("rbind", parseArgs(input)))
inputL <- as.list(as.character(inputDF$V2))
names(inputL) <- inputDF$V1

# Original Script

library(rjags)
library(data.table)

dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=fread("Lawless_area.txt",header=FALSE)
    times=fread("Lawless_time.txt",header=FALSE)
    data=fread("Lawless_data.txt",header=FALSE) #3rd column (Identifier) => strain_parentcolony 
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=fread("Levy_area.txt",header=FALSE)
    times=fread("Levy_time.txt",header=FALSE)
    data=fread("Levy_data.txt",header=FALSE) #3rd column (Identifier) => replicate
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    area=fread("Ziv_area_filtered1.txt",header=FALSE)
    times=fread("Ziv_times_filtered1.txt",header=FALSE)
    data=fread("Ziv_data_filtered1.txt",header=FALSE) #3rd column (Identifier) => colony
    return(list("area"=area,"data"=data,"times"=times))
  }
  else {print("Not a valid dataset")}
}

datsetname=inputL$arg1

# Choosing a data set 
x=dataset(datsetname)
area=x$area
times=x$times
data=x$data

#Parameters for prior distributions
intercept=max(area$V1) # this is the maximum value for x0
kval1=max(area$V23,na.rm=TRUE)
kval2=min(area$V23,na.rm=TRUE)

modelstring="
model{
  for (i in 1:n) {
    mu[i]<-(K*x0*exp(r*t[i]))/(K+x0*(exp(r*t[i])-1))
    x[i]~dnorm(mu[i],tau)
  }
  for (j in 1:n) {
    pmu[j]<-(K*x0*exp(r*t[j]))/(K+x0*(exp(r*t[j])-1))
    px[j]~dnorm(pmu[j],tau)
  }
  K~dunif(0,2000000)
  r~dunif(0,2)
  x0~dunif(0,2000)
  tau~dunif(0,1000)  
}
"

modelstring2="
model{
  for (i in 1:n){
    mu[i]<-x0*exp(r*t[i])
    x[i]~dnorm(mu[i],tau)
  }
  for (j in 1:n) {
    pmu[j]<-x0*exp(r*t[j])
    px[j]~dnorm(pmu[j],tau)
  }
  r~dunif(0,2)
  x0~dunif(0,2000)
  tau~dunif(0,1000)  
}
"

pdf(height = 16, width = 16, file = paste(datsetname,"_Logistic_Model_Output_",inputL$arg2,".pdf",sep=""))
#pdf(height = 16, width = 16, file = paste(datsetname,"_Exponential_Model_Output.pdf",sep=""))

x0_total=c()
r_total=c()
k_total=c()

growth_curves=strsplit(inputL$arg2,":")
growth_curves=as.numeric(noquote(growth_curves[[1]][1])):as.numeric(noquote(growth_curves[[1]][2]))

for (i in growth_curves) { #1:dim(area)[1]
  print(i)
  vals=which(!is.na(area[1,]))
  N=length(vals)
  prednames=sprintf("pmu[%i]",1:N)
  dat=list('x'=as.numeric(area[i,])[vals], 't'=as.numeric(times[i,])[vals], 'n'=N)
  
  jags<-jags.model(textConnection(modelstring),
                   data=dat, #names must be those in the JAGS model specification
                   n.chains=4) #how many parallele chains to run
                   #n.adapt=1000) #how many samples to throw away as part of the adaptive sampling period for each chain
  update(jags,1000000) #Burn-in period
  
  samples=coda.samples(model=jags,variable.names=c('r','x0','tau','K', prednames),n.iter=1000000,thin=1000)
  subset=samples[,c("r","K","x0","tau")]
  
  # samples=coda.samples(model=jags,variable.names=c('r','x0','tau', prednames),n.iter=1000000,thin=1000)
  # subset=samples[,c("r","x0","tau")]
  
  plot(subset)
  #plot(samples)
  #print(summary(samples))

  op=par(mfrow=c(2,2))
  #op=layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  
  #Plot the Growth Curve 
  plot(as.numeric(times[i,])[vals],as.numeric(area[i,])[vals],type='l',xlab="Time",ylab="Area",lty=3, main=paste("Growth Curve",i))
  points(as.numeric(times[i,])[vals],as.numeric(area[i,])[vals],type='p',pch=16)
  #Overlaying the Posterior Predictive 
  preds=samples[[1]][,prednames]
  lines(as.numeric(times[i,])[vals],apply(preds,2,mean),type="l",col="red",ylim=c(0,1000),xlab="Time",ylab="Population Size",main="Posterior predictive")
  points(as.numeric(times[i,])[vals],apply(preds,2,quantile,0.05),type="l",lty=3,col="red")
  points(as.numeric(times[i,])[vals],apply(preds,2,quantile,0.95),type="l",lty=3,col="red")
  
  #Posterior Distribution for Growth Rate, r
  m_r=summary(samples)[[1]]['r',1]
  print(paste("Mean of r is", m_r))
  curve(dunif(x,0,2),from=0,to=2,main=paste("Growth Rate, r; Mean:",signif(m_r,3)),ylim=c(0,max(density(samples[[1]][,"r"])$y)))
  points(density(samples[[1]][,"r"]),type="l",col="blue")
  abline(v=m_r,col="black",lty=3)
  
  #Posterior Distribution for Intercept, x0
  m_x0=summary(samples)[[1]]['x0',1]
  print(paste("Mean of x0 is", m_x0))
  curve(dunif(x,0,2000),from=0,to=2000,main=paste("Intercept, x0; Mean:",signif(m_x0,3)),ylim=c(0,max(density(samples[[1]][,"x0"])$y)))
  points(density(samples[[1]][,"x0"]),type="l",col="blue")
  abline(v=m_x0,col="black",lty=3)
  
  #Posterior Distribution for Carrying Capacity, K
  m_k=summary(samples)[[1]]['K',1]
  print(paste("Mean of k is", m_k))
  curve(dunif(x,0,2000000),from=0,to=2000000,main=paste("Carrying Capacity, K; Mean:",signif(m_k,3)),ylim=c(0,max(density(samples[[1]][,"K"])$y)))
  points(density(samples[[1]][,"K"]),type="l",col="blue")
  abline(v=m_k,col="black",lty=3)
  
  par(op)

  x0_total=c(x0_total,m_x0)
  r_total=c(r_total,m_r)
  k_total=k_total=c(k_total,m_k)

}

dev.off()

Bayes_parameters=data.frame("Intercept"=x0_total,"Rate"=r_total,"CarryingCapacity"=k_total)
filename=paste(datsetname,"_Bayes_parameters_",inputL$arg2,".txt",sep="")
write.table(Bayes_parameters,filename,col.names=TRUE,row.names=FALSE)