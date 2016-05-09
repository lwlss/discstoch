# Fitting three different variations of the exponential model to the data to see which one works best: 

library(data.table)

# Choosing a data set to work with
dataset<-function(x){
  if (x == "Lawless"){
    # DataSet1: Lawless
    area=fread("Lawless_area.txt",header=FALSE)
    times=fread("Lawless_time.txt",header=FALSE)
    data=fread("Lawless_data.txt",header=TRUE) #3rd column (Identifier) => strain_parentcolony 
    residuals=fread("Lawless_residuals.txt",header=FALSE)
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    area=fread("Levy_area.txt",header=FALSE)
    times=fread("Levy_time.txt",header=FALSE)
    data=fread("Levy_data.txt",header=TRUE) #3rd column (Identifier) => replicate
    residuals=fread("Levy_residuals.txt",header=FALSE)
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    area=fread("Ziv_area.txt",header=FALSE)
    times=fread("Ziv_time.txt",header=FALSE)
    data=fread("Ziv_data.txt",header=TRUE) #3rd column (Identifier) => colony
    residuals=fread("Ziv_residuals.txt",header=FALSE)
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else {print("Not a valid dataset")}
}

datsetname="Lawless"
x=dataset(datsetname)
area=as.matrix(x$area)
times=as.matrix(x$times)
data=x$data
names(data)=c("genotype","clonalcolony","identifier")
area[area==1]=NA #Lawless only

# Choosing growth curves which must have at least 10 or more consecutive datapoints (Lawless)
selected_indices=c()
for (i in 1:dim(area)[1]){
  check_indices=which(area[i,]>0) #this does not take the NA values into account 
  testing=rle(diff(check_indices))
  if (any(testing$lengths>=10 & testing$values==1) == TRUE){
    selected_indices=c(selected_indices, i)
  } #checking whether there are at least 3 consecutive indices 
}

area=area[selected_indices,]
times=times[selected_indices,]
data=data[selected_indices,]

area=area[-4030,]
area=area[-4743,]
times=times[-4030,]
times=times[-4743,]
data=data[-4030,]
data=data[-4743,]

write.table(area,paste(datsetname,"_area_filtered1.txt",sep=""),col.names=FALSE,row.names=FALSE)
write.table(times,paste(datsetname,"_times_filtered1.txt",sep=""),col.names=FALSE,row.names=FALSE)
write.table(data,paste(datsetname,"_data_filtered1.txt",sep=""),col.names=TRUE,row.names=FALSE)

# ################################### Model 1 ###########################################################
# # The Current Exponential Model
# AList=c()
# rList=c()
# residuals=matrix(NA,nrow=dim(area)[1],ncol=dim(area)[2]) #residuals maybe of different lengths; thus NA for when no data point is available
# for (i in 1:dim(area)[1]){
#   if(length(unique(area[i,]))>1){
#     # Getting rid of NA values
#     indices=as.numeric(which(area[i,]>0))
#     i_area=as.numeric(area[i,indices])
#     i_times=as.numeric(times[i,indices])
#     # Fit exponential model
#     result=lm(y~x,data=data.frame(x=i_times,y=log(i_area)))
#     resid=residuals(result)
#     A=exp(result$coefficients[[1]]) #converting back from logs
#     r=(result$coefficients[[2]]) #this still includes negative growth rates which have no biological meaning
#   }else{A=1;r=0;residuals=rep(NA,dim(area)[2])}
#   AList=c(AList,A)
#   rList=c(rList,r)
#   residuals[i,indices]=resid
# }
# 
# data$rate=rList
# data$intercept=AList
# 
# write.table(cbind(data,residuals),file=paste(datsetname,"_data_model1.txt",sep=""))


################################### Model 2 ###########################################################

# Fitting the logarithmic exponential model using nls
AList=c()
rList=c()
residuals=matrix(NA,nrow=dim(area)[1],ncol=dim(area)[2]) #residuals maybe of different lengths; thus NA for when no data point is available
for (i in 1:dim(area)[1]){
  if(length(unique(area[i,]))>1){
    inocguess=min(area[i,],na.rm=TRUE) #minimum data point;
    # Getting rid of NA values
    indices=as.numeric(which(area[i,]>0))
    i_area=as.numeric(area[i,indices])
    i_times=as.numeric(times[i,indices])
    # Fit exponential model
    result=nls(y~A+r*x,data=data.frame(x=i_times,y=log(i_area)),start=list(A=log(inocguess),r=0.1))
    resid=residuals(result)
    res=result$m$getPars()
    A=exp(as.numeric(res[1])) #converting back from logs
    r=(as.numeric(res[2])) #this still includes negative growth rates which have no biological meaning
  }else{A=1;r=0;residuals=rep(NA,dim(area)[2])}
  AList=c(AList,A)
  rList=c(rList,r)
  residuals[i,indices]=resid
}

data$rate=rList
data$intercept=AList

write.table(cbind(data,residuals),file=paste(datsetname,"_data_model2.txt",sep=""),col.names=TRUE,row.names=FALSE)
# 
# 
# ################################### Model 3 ###########################################################
# 
# 
# # Fitting the exponential model using nls
# 
# #Lawless Error in nls(y ~ A * exp(r * x), data = data.frame(x = i_times, y = log(i_area)),  :
# #step factor 0.000488281 reduced below 'minFactor' of 0.000976562
# 
# pdf(height = 16, width = 16, file = "Exponential_Model_Fit.pdf")
# par(mfrow=c(4,4))
# 
# AList=c()
# rList=c()
# residuals=matrix(NA,nrow=dim(area)[1],ncol=dim(area)[2]) #residuals maybe of different lengths; thus NA for when no data point is available
# for (i in 1:dim(area)[1]){
#   print(i)
#   if(length(unique(area[i,]))>1){
#     inocguess=min(area[i,],na.rm=TRUE) #minimum data point;
#     # Getting rid of NA values
#     indices=as.numeric(which(area[i,]>0))
#     i_area=as.numeric(area[i,indices])
#     i_times=as.numeric(times[i,indices])
#     # Fit exponential model
#     result=nls(y~A*exp(r*x),data=data.frame(x=i_times,y=i_area),start=list(A=inocguess,r=0.1),model=TRUE)
#     resid=residuals(result)
#     res=result$m$getPars()
#     A=(as.numeric(res[1]))
#     r=(as.numeric(res[2])) #this still includes negative growth rates which have no biological meaning
#     # Producing the plot 
#     plot(predict(result),type='l',main=paste("Growth Curve:", i),ylab="Area",xlab="Time",cex.main=1.4,cex.lab=1.4)
#     lines(area[i,],type='p')
#   }else{A=1;r=0;residuals=rep(NA,dim(area)[2])}
#   AList=c(AList,A)
#   rList=c(rList,r)
#   residuals[i,indices]=resid
# }
# 
# dev.off()
# 
# 
# data$rate=rList
# data$intercept=AList
# 
# write.table(cbind(data,residuals),file=paste(datsetname,"_data_model3.txt",sep=""),col.names=TRUE,row.names=FALSE)
