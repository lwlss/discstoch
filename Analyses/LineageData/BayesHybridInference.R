library(smfsb)
library(detstocgrowth)

dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=fread("~/BayesianInference/Lawless_area_shortTC.txt",header=FALSE)
    times=fread("~/BayesianInference/Lawless_time_shortTC.txt",header=FALSE)
    data=fread("~/BayesianInference/Lawless_data_shortTC.txt",header=FALSE) #3rd column (Identifier) => strain_parentcolony
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=fread("~/BayesianInference/Levy_area_filtered1.txt",header=FALSE)
    times=fread("~/BayesianInference/Levy_times_filtered1.txt",header=FALSE)
    data=fread("~/BayesianInference/Levy_data_filtered1.txt",header=FALSE) #3rd column (Identifier) => replicate
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    area=fread("~/BayesianInference/Ziv_area_filtered1.txt",header=FALSE)
    times=fread("~/BayesianInference/Ziv_times_filtered1.txt",header=FALSE)
    data=fread("~/BayesianInference/Ziv_data_filtered1.txt",header=FALSE) #3rd column (Identifier) => colony
    return(list("area"=area,"data"=data,"times"=times))
  }
  else {print("Not a valid dataset")}
}

# Choosing a data set
datsetname="Lawless"

# Extracting data for the growth curve(s) on which to do inference
gc=252
x=dataset(datsetname)
area=x$area[gc,]
times=x$times[gc,]
data=x$data[gc,]
# Apply calibration and format to feed to the Hybrid Model
cells=as.numeric(area)/area[,1]
data=data.frame(c=cells,t=as.numeric(times))
rownames(data)=data$t
data$t=NULL
# Plot the Growth Curve
plot(rownames(data),data$c,xlab="Time (h)",ylab="Cell number",cex.lab=1.5,main="Growth Curve",
     type='p',pch=16,col="red")

# Specify parameters for Inference
guessK=100000
guessr=0.3
params=c(K=guessK,r=guessr)
pmin=c(K=0.5*guessK,r=0.01)
pmax=c(K=2*guessK,r=1)

# Inference using the Hybrid Model
mcmc_chain=BayesHybrid(params,pmin,pmax,1000000,0.02,1000,20,10,data)
mcmcSummary(mcmc_chain,show=TRUE,plot=TRUE)

# Plotting the MCMC Output
PlotPosteriorProb(pmin,pmax,mcmc_chain)

# Plotting the Posterior Predictives
PlotPosteriorPredictive(data,parmat,10,guessK)

# Write first MCMC chain to file
write.table(mcmc_chain,"GC252HybBayesLog.txt",col.names=TRUE,row.names=FALSE)
