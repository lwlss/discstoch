#ADD-HOC ANALYSES 

library(nortest)
library(moments)
library(data.table)

############################################################################################################

# Analysing the lag phase at the population level
#x=fread("PopSimStats.txt",header=TRUE)

############################################################################################################

# Analysing the residuals at the single lineage level
datset="Levy"
readin=paste(datset,"_residuals.txt",sep="")

residuals=fread(readin)

# Analysing the mean residual fit for the exponential model
meanRes=rowSums(residuals, na.rm=TRUE)/dim(residuals)[2]
# op=par(mfrow=c(1,1))
# total_mean=mean(meanRes)
# plot(meanRes, main=paste("Mean Residuals Values for Each Lineage, N = ",toString(dim(residuals)[1]),sep=""))
# par(op)
print(mean(meanRes))

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

normvals=testing_normality(residuals)
print(dim(residuals)[1])
print(length(which(normvals$Lilliefors.H0==TRUE)))
print(length(which(normvals$ShapiroWilk.H0==TRUE)))

#Storing the output in a file
op_residuals=cbind(data.frame(normvals),residuals)
writeout=paste(datset,"_residuals_normtests.txt",sep="")
write.table(op_residuals,writeout,col.names=TRUE,row.names=FALSE)
x=fread(writeout,header=TRUE)

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

