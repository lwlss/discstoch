#Adding a calculated growth rate value to each of the data set 
#In order to speed up later calculations 

library(data.table)

#################Functions################ (prevous functions re-used)

dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=fread("Lawless_area.txt",header=FALSE)
    times=fread("Lawless_time.txt",header=FALSE)
    data=fread("Lawless_data.txt",header=TRUE) #3rd column (Identifier) => strain_parentcolony 
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=fread("Levy_area.txt",header=FALSE)
    times=fread("Levy_time.txt",header=FALSE)
    data=fread("Levy_data.txt",header=TRUE) #3rd column (Identifier) => replicate
    return(list("area"=area,"data"=data,"times"=times))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    area=fread("Ziv_area.txt",header=FALSE)
    times=fread("Ziv_time.txt",header=FALSE)
    data=fread("Ziv_data.txt",header=TRUE) #3rd column (Identifier) => colony
    return(list("area"=area,"data"=data,"times"=times))
  }
  else {print("Not a valid dataset")}
}

#Calculating the growth rate (LS for exponential model)
LM_growthrate<-function(ai,ti){
  fit<-lm(log(ai)~ti)
  res=resid(fit)
  rate=fit$coefficient[[2]]
  intercept=fit$coefficients[[1]]
  return(list("rate"=rate,"res"=res,"int"=intercept))
}

#######################Main###################
# Choosing a data set 
x=dataset("Ziv") #adjust the data set accordingly
area=x$area
times=x$times
data=x$data
area[area==0]=NA

# #Obtaining the growth rate
# k=c()
# residuals=matrix(NA,nrow=dim(area)[1],ncol=dim(area)[2])
# for (i in 1:dim(area)[1]){
#   rate=LM_growthrate(as.numeric(area[i,]),as.numeric(times[i,]))$rate
#   res=LM_growthrate(as.numeric(area[i,]),as.numeric(times[i,]))$res
#   residuals[i,as.numeric(names(res))]=res
# }

#Obtaining the intercept 
intercepts=c()
for (i in 1:dim(area)[1]){
  inter=LM_growthrate(as.numeric(area[i,]),as.numeric(times[i,]))$int
  intercepts=c(intercepts,inter)
}


# #Storing the growth rate 
# data$rate=k

#Storing the intercept
data$intercept=intercepts


#Writing to the file 
write.table(data,"Ziv_test1.txt",col.names=TRUE,row.names=FALSE)
#read.table("Lawless_test1.txt,header=TRUE")
#write.table(residuals,"Levy_residuals.txt",col.names=FALSE,row.names=FALSE)

