# Bayesian Deterministic Parameter Inference

library(detstocgrowth)
library(data.table)

dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=fread("~/Lawless_area_shortTC.txt",header=FALSE)
    times=fread("~/Lawless_time_shortTC.txt",header=FALSE)
    data=fread("~/Lawless_data_shortTC.txt",header=FALSE) #3rd column (Identifier) => strain_parentcolony
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=fread("~/Levy_area_filtered1.txt",header=FALSE)
    times=fread("~/Levy_times_filtered1.txt",header=FALSE)
    data=fread("~/Levy_data_filtered1.txt",header=FALSE) #3rd column (Identifier) => replicate
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    area=fread("~/Ziv_area_filtered1.txt",header=FALSE)
    times=fread("~/Ziv_times_filtered1.txt",header=FALSE)
    data=fread("~/Ziv_data_filtered1.txt",header=FALSE) #3rd column (Identifier) => colony
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

# Write output to PDF file
pdf(height = 16, width = 16, file = paste(datsetname,"_Exponential_Model_Output.pdf",sep=""))
bayesinf=BayesDet(area,times,2,1000000,1000,2)
dev.off()

# Write first MCMC chain to file
write.table(bayesinf$Samples[[1]],"GC252DetBayesExp.txt",col.names=TRUE,row.names=FALSE)
