#ADD-HOC ANALYSES 

library(nortest)
library(moments)
library(data.table)
library(car)
library(lmtest)

############################################################################################################

# Analysing the lag phase at the population level
#x=fread("PopSimStats.txt",header=TRUE)

############################################################################################################

# Analysing the residuals at the single lineage level

####Functions####
dataset<-function(x,ou.rm=FALSE){ #can choose to take area where outliers have been removed
  if (x == "Lawless"){
    # DataSet1: Lawless
    if(ou.rm==FALSE){area=fread("Lawless_area.txt",header=FALSE)}else{area=fread("Lawless_area_or.txt",header=FALSE)}
    times=fread("Lawless_time.txt",header=FALSE)
    data=fread("Lawless_data.txt",header=TRUE) #3rd column (Identifier) => strain_parentcolony 
    residuals=fread("Lawless_residuals.txt",header=FALSE)
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Levy"){
    # DataSet2: Levy
    if(ou.rm==FALSE){area=fread("Levy_area.txt",header=FALSE)}else{area=fread("Levy_area_or.txt",header=FALSE)}
    times=fread("Levy_time.txt",header=FALSE)
    data=fread("Levy_data.txt",header=TRUE) #3rd column (Identifier) => replicate
    residuals=fread("Levy_residuals.txt",header=FALSE)
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else if (x == "Ziv"){
    # DataSet3: Ziv
    if(ou.rm==FALSE){area=fread("Ziv_area.txt",header=FALSE)}else{area=fread("Ziv_area_or.txt",header=FALSE)}
    times=fread("Ziv_time.txt",header=FALSE)
    data=fread("Ziv_data.txt",header=TRUE) #3rd column (Identifier) => colony
    residuals=fread("Ziv_residuals.txt",header=FALSE)
    return(list("area"=area,"data"=data,"times"=times,"residuals"=residuals))
  }
  else {print("Not a valid dataset")}
}
subset_strain<-function(d,a,t,strain){ 
  s_name=toString(strain)
  indices=which(d$genotype == strain)
  s_area=a[indices,]; s_times=t[indices,]; s_data=d[indices,]
  return(list("area"=s_area,"times"=s_times,"data"=s_data, "name"=s_name, "indices"=indices))
}
subset_colony<-function(d,a,t,colony){
  c_name=toString(colony)
  indices=which(d$clonalcolony == colony)
  c_area=a[indices,]; c_times=t[indices,]; c_data=d[indices,]
  return(list("area"=c_area,"times"=c_times,"data"=c_data, "name"=c_name))
}
subset_identifier<-function(d,a,t,identifier){
  i_name=toString(identifier)
  indices=which(d$identifier == identifier)
  i_area=a[indices,]; i_times=t[indices,]; i_data=d[indices,]
  return(list("area"=i_area,"times"=i_times,"data"=i_data, "name"=i_name))
}
testing_normality<-function(r){
  tskew=c()
  tktosis=c()
  tpval1=c()
  tpval2=c()
  tTF1=c()
  tTF2=c()
  for(i in 1:dim(r)[1]){
    print(i)
    ri=as.numeric(r[i,])
    if(sum(ri,na.rm=TRUE)==0){
      skew=NA
      ktosis=NA
      pval1=NA
      TF1=TRUE #don't want to discard a perfect fit; although this will always be a straight line (all 1s or all 0s)
      pval2=NA
      TF2=TRUE 
      tskew=c(tskew,skew)
      tktosis=c(tktosis,ktosis)
      tpval1=c(tpval1,pval1)
      tpval2=c(tpval2,pval2)
      tTF1=c(tTF1,TF1)
      tTF2=c(tTF2,TF2)
    }
    else if(sum(!is.na(residuals[3,]))<=4){
      skew=NA
      ktosis=NA
      pval1=NA
      TF1=FALSE #lineages of 4 or less data points will be discarded! 
      pval2=NA
      TF2=FALSE 
      tskew=c(tskew,skew)
      tktosis=c(tktosis,ktosis)
      tpval1=c(tpval1,pval1)
      tpval2=c(tpval2,pval2)
      tTF1=c(tTF1,TF1)
      tTF2=c(tTF2,TF2)
    }
    else{
      skew=skewness(ri,na.rm=TRUE)
      ktosis=kurtosis(ri,na.rm=TRUE)
      test1=lillie.test(ri) #H0: data are normally distributed; p>0.1: we do not reject H0
      test2=shapiro.test(ri)
      pval1=test1$p.value
      if (pval1<0.1){
        TF1=FALSE
      }else{TF1=TRUE}
      pval2=test2$p.value
      if (pval2<0.1){
        TF2=FALSE
      }else{TF2=TRUE}
      tskew=c(tskew,skew)
      tktosis=c(tktosis,ktosis)
      tpval1=c(tpval1,pval1)
      tpval2=c(tpval2,pval2)
      tTF1=c(tTF1,TF1) #p>0.1
      tTF2=c(tTF2,TF2)
    }
  }
  return(list("skewness"=tskew,"kurtosis"=tktosis,"Lilliefors.p"=tpval1, "Lilliefors.H0"=tTF1, "ShapiroWilk.p"=tpval2, "ShapiroWilk.H0"=tTF2))
}
LevTest<-function(pick){
  #Formatting the data
  area_list=c(t(pick$area[1:dim(pick$area)[1],]))
  categor=matrix(0,nrow=dim(pick$area)[1],ncol=dim(pick$area)[2])
  for (i in 1:dim(pick$area)[1]){
    categor[i,]=c(rep(i,dim(pick$area)[2]))
  }
  categor_list=factor(c(t(categor)))
  #Testing for equality of variance in "pick" population
  lev_data=data.frame("area_x"=log(area_list),"lineage"=categor_list)
  ltest=leveneTest(lm(area_x~lineage,data=lev_data),na.rm=TRUE)
  return(ltest)
}

######Data#####
datsetname="Lawless"
x=dataset(datsetname,ou.rm=TRUE)
area=x$area
times=x$times
data=x$data
residuals=x$residuals
area[area==0]=NA
identifier_names=unique(data$identifier)

#Variances among strains (even those stemming from the parent paretn colony) are not equal... Heterogeneity in growth rate! 
#Can use this test to justify why we are looking at single lineages, maybe... results are somewhat suspiciously close to zero
LT_p=c()
for(i in 1:length(identifier_names)){
  pickid=identifier_names[i] #choose strain here!
  identifier=subset_identifier(data,or_area,times,pickid)
  LT=LevTest(identifier)
  LT_p=c(LT_p,LT$`Pr(>F)`)
}

###############StoredData#################

# # Analysing the mean residual fit for the exponential model
# meanRes=rowSums(residuals, na.rm=TRUE)/dim(residuals)[2]
# # op=par(mfrow=c(1,1))
# # total_mean=mean(meanRes)
# # plot(meanRes, main=paste("Mean Residuals Values for Each Lineage, N = ",toString(dim(residuals)[1]),sep=""))
# # par(op)
# print(mean(meanRes))


# #Removing outliers with a Cook's distance >1 
# #Calculating the Breusch Pagan values based on that area 
# 
# or_area=area #area with outliers removed 
# BreuschPagan=c()
# for (i in 1:dim(area)[1]){
#   model=lm(log(as.numeric(area[i,]))~as.numeric(times[i,]))
#   if (sum(as.numeric(as.numeric(cooks.distance(model)))>1,na.rm=TRUE) >= 1){
#     or_area[i,which(as.numeric(cooks.distance(model))>1)] = NA #eliminate all data point with CD>1
#     print(i); print(which(as.numeric(cooks.distance(model))>1))
#   }
#   BreuschPagan=c(BreuschPagan,as.numeric(bptest(model)$p.value)) #NOTE, this return NaN for when all values are the same 
# }
# 
# length(which(BreuschPagan<0.05)) # quite a lot but less than when testing for normality.... 
# 
# #Storing the outputs in a file 
# overwrite=paste(datsetname,"_residuals_normtests.txt",sep="")
# x=fread(overwrite,header=TRUE)
# write.table(cbind(BreuschPagan,x),overwrite,col.names=TRUE,row.names=FALSE)
# write.table(or_area,paste(datsetname,"_area_or.txt",sep=""),col.names=FALSE,row.names=FALSE)


# normvals=testing_normality(residuals)
# print(dim(residuals)[1])
# print(length(which(normvals$Lilliefors.H0==TRUE)))
# print(length(which(normvals$ShapiroWilk.H0==TRUE)))
# 
# #Storing the output in a file
# op_residuals=cbind(data.frame(normvals),residuals)
# writeout=paste(datset,"_residuals_normtests.txt",sep="")
# write.table(op_residuals,writeout,col.names=TRUE,row.names=FALSE)
# x=fread(writeout,header=TRUE)


############################################################################################################

# #Analysing  the intercept values of the exponential model 
# datset="Levy"
# readin=paste(datset,"_test1.txt",sep="")
# dat=fread(readin,header=TRUE)
# int=dat$intercept
# 
# hist_cells<-function(dat,interval){
#   lo=trunc(min(dat)*10)/10-0.1 #rounding down
#   hi=trunc(max(dat)*10)/10+0.1 #rounding up 
#   cells=seq(lo,hi,interval)
#   return(cells)
# }
# 
# cells=hist_cells(int,0.1)
# op=par(mfrow=c(1,1))
# histo=hist(int, breaks=cells, xlab="log(Area)", main=paste(datset, "; Intercept Values for the Exponential Model", sep=""),cex.lab=1.4, cex.main=1.2,col="lightblue")
# valguess=(histo$breaks[which.max(histo$counts)]+histo$breaks[which.max(histo$counts)+1])/2
# abline(v=valguess,col="red",cex=1.5)
# abline(v=median(int),col="green")
# abline(v=mean(int),col="darkorange")
# legend("topleft",legend=c("LME guess","Median","Mean"), col=c("red","green","darkorange"),pch=20)
# par(op)
# 
# valguess=(histo$breaks[which.max(histo$counts)]+histo$breaks[which.max(histo$counts)+1])/2
# print(paste("Estimated area of one cell is", valguess))
# print(paste("Median starting area is",median(int)))
# print(paste("Mean starting area is",mean(int)))

