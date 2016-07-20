# Analysing Noise in Growth Curves 

library(data.table)

dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=fread("~/BayesianInference/Lawless_area_shortTC.txt",header=FALSE)
    times=fread("~/BayesianInference/Lawless_time_shortTC.txt",header=FALSE)
    data=fread("~/BayesianInference/Lawless_data_shortTC.txt",header=FALSE) #3rd column (Identifier) => strain_parentcolony
    names(data)=c("genotype","clonalcolony","identifier","blobnumber")
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=fread("~/BayesianInference/Levy_area_filtered1.txt",header=FALSE)
    times=fread("~/BayesianInference/Levy_times_filtered1.txt",header=FALSE)
    data=fread("~/BayesianInference/Levy_data_filtered1.txt",header=TRUE) #3rd column (Identifier) => replicate
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    area=fread("~/BayesianInference/Ziv_area_filtered1.txt",header=FALSE)
    times=fread("~/BayesianInference/Ziv_times_filtered1.txt",header=FALSE)
    data=fread("~/BayesianInference/Ziv_data_filtered1.txt",header=TRUE) #3rd column (Identifier) => colony
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
area=as.matrix(area)

noise=matrix(0,nrow=dim(area)[1],ncol=dim(area)[2]-1)

for (i in 1:dim(area)[1]){
  noise[i,]=abs(as.numeric(diff(area[i,])))
}

meanNoise=colSums(noise)/dim(noise)[1]

print(meanNoise)

strain_names=unique(data$genotype)
his3=detstocgrowth::subset_strain(data,area,times,strain_names[1])
htz1=detstocgrowth::subset_strain(data,area,times,strain_names[2])

meanhis3=colSums(noise[his3$indices,])/length(his3$indices)
meanhtz1=colSums(noise[htz1$indices,])/length(htz1$indices)

print(meanhis3)
print(meanhtz1)