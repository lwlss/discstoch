# Bayesian Parameter Inference
# Set working directory to source file location 

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
    area=fread("Ziv_area.txt",header=FALSE)
    times=fread("Ziv_time.txt",header=FALSE)
    data=fread("Ziv_data.txt",header=FALSE) #3rd column (Identifier) => colony
    return(list("area"=area,"data"=data,"times"=times))
  }
  else {print("Not a valid dataset")}
}

datsetname="Lawless"

# Choosing a data set 
x=dataset(datsetname)
area=x$area
times=x$times
data=x$data
area[area==0]=NA

N=dim(area)[2]

pdf(height = 16, width = 16, file = "Linear_Regression_Output_Plots_1.pdf")
a_total=c()
b_total=c()
for (i in 1:dim(area)[1]){ #64:68 for example works fine
  print(i)
  jags<-jags.model('Linear_Regression_1.bug',
                   data=list('y'=log(as.numeric(area[i,])), 'x'=as.numeric(times[i,]), 'N'=N), #names must be those in the JAGS model specification
                   n.chains=4, #how many parallele chains to run
                   n.adapt=1000) #how many samples to throw away as part of the adaptive sampling period for each chain
  update(jags,100) #Burn-in period
  op1=jags.samples(jags,c('a','b'),1000)
  samples=coda.samples(jags,c('a','b'),1000)
  plot(samples)
  #plot(samples[,'b']) #just looking at the estimated growth rate 
  a_total=c(a_total,mean(op1$a))
  b_total=c(b_total,mean(op1$b))
}

dev.off()

Bayes_parameters=data.frame("Intercept"=a_total,"Rate"=b_total)
filename=paste(datsetname,"_Bayes_parameters.txt",sep="")
write.table(Bayes_parameters,filename,col.names=TRUE,row.names=FALSE)


