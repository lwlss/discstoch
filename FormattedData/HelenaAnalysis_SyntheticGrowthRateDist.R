popsimdat=function(k,N,t,it){
  PopData=matrix(0,nrow=it,ncol=length(t))
  total_indices=c()
  total_simdata=list()
  total_rates=c()
  for (i in 1:it){
    indices=sample(1:length(k),N,replace=TRUE)
    rates=k[indices]
    simdata=matrix(0,nrow=length(rates),ncol=length(t))
    for (j in 1:length(rates)){
      simdata[j,]=simExponential(rates[j],t)
    }
    total_simdata[[i]]=simdata
    PopData[i,]=colSums(simdata)
    total_indices=c(total_indices,indices)
    total_rates=c(total_rates,rates)
  }
  return(list("PopData"=PopData,"SimData"=total_simdata,
              "indices"=total_indices,"SimRates"=total_rates))
}

N=100
iter=10000
time=seq(0,35,0.5)

plotpopsimdat=function(popdata,r,time,iterations){
  pop_rates=c()
  for (i in 1:dim(popdata)[1]){
    k=detstocgrowth::LM_growthrate(as.numeric(popdata[i,15:35]),time)$rate
    pop_rates=c(pop_rates,k)
  }
  plot(1,type='n', xlim=range(time), ylim=range(popdata),xlab="Time (h)",
       ylab="No. of Cells",log='y',cex.lab=1.4)
  colours=rainbow(20)
  # colours=colours[-1]
  # colours=colours[-1]
  gr_range=seq(0.05,1,0.05)
  #gr_range=seq(0.05,0.4,0.05)
  # bp_pr=c()
  # bp_loc=c()
  for (i in 1:iterations){
    id=which(abs(gr_range-max(poprate[(i*100-99):(i*100-100+N)]))==
               min(abs(gr_range-max(poprate[(i*100-99):(i*100-100+N)]))))
    lines(time,popdata[i,],col=colours[id],lwd=2)
    # if(max(poprate[(i*100-99):(i*100-100+N)])>0.3){ 
    #   op_bcp=bcp(log(popdata[i,]),time)
    #   max_prob=max(op_bcp$posterior.prob,na.rm=TRUE)
    #   breakpoint_location=which(op_bcp$posterior.prob==max_prob)
    #   bp_pr=c(bp_pr,max_prob)
    #   bp_loc=c(bp_loc,breakpoint_location)
    # }
  }
  # print(mean(bp_pr))
  # print(min(bp_pr))
  # print(time[mean(bp_loc)])
  # legend("topleft",title="Max. Growth Rate",legend=rev(gr_range)[1:3],
  #        col=rev(colours)[1:3],lty=rep(1,3),lwd=rep(3,3))
  
  h=hist(r,breaks=seq(0,ceiling(max(r*10))/10,0.01),plot=FALSE)
  hcol=colours[sapply(h$mids, function(x) which(abs(gr_range-x)==min(abs(gr_range-x)))[1])]
  plot(h,col=hcol,cex.lab=1.4,main="")
  abline(v=mean(poprate),col="black",lwd=3)
  abline(v=mean(pop_rates),col="black",lwd=3,lty=3)
  # legend("topright",legend=c("True Mean","Pop'n Mean"),
  #        col=c("black","black"),lty=c(1,2),lwd=c(3,3))
}

op=par(mfrow=c(5,2))
norm=rnorm(1000,0.2,0.04)
norm[which(norm<0)]=0
# hist(norm,breaks=seq(0,ceiling(max(norm*10))/10,0.01),
#      cex.lab=1.4,xlab="Growth Rate (1/h)",main="")
r=norm
popsim=popsimdat(r,N,time,iter)
popdata=popsim$PopData
poprate=popsim$SimRates
plotpopsimdat(popdata,r,time,iter)

widenorm=rnorm(1000,0.2,0.1)
widenorm[which(widenorm<0)]=0
widenorm[which(widenorm>0.4)]=rnorm(length(which(widenorm>0.4)),0.2,0.001)
# hist(widenorm,breaks=seq(0,ceiling(max(widenorm*10))/10,0.01),
#      cex.lab=1.4,xlab="Growth Rate (1/h)",main="")
r=widenorm
popsim=popsimdat(r,N,time,iter)
popdata=popsim$PopData
poprate=popsim$SimRates
plotpopsimdat(popdata,r,time,iter)

longtail=norm[-sample(length(norm),10)]
longtail=c(longtail,rnorm(10,0.75,0.1))
# hist(longtail,breaks=seq(0,ceiling(max(longtail*10))/10,0.01),
#      cex.lab=1.4,xlab="Growth Rate (1/h)",main="")
r=longtail
popsim=popsimdat(r,N,time,iter)
popdata=popsim$PopData
poprate=popsim$SimRates
plotpopsimdat(popdata,r,time,iter)

bimode=rnorm(600,0.2,0.05)
bimode=c(bimode,rnorm(400,0.4,0.075))
bimode[which(bimode<0)]=0
# hist(bimode,breaks=seq(0,ceiling(max(bimode*10))/10,0.01),
#      cex.lab=1.4,xlab="Growth Rate (1/h)",main="")
r=bimode
popsim=popsimdat(r,N,time,iter)
popdata=popsim$PopData
poprate=popsim$SimRates
plotpopsimdat(popdata,r,time,iter)

trimode=bimode[-sample(length(bimode),100)]
trimode=c(trimode,rnorm(100,0,0.1))
trimode[which(trimode<0)]=0
# hist(trimode,breaks=seq(0,ceiling(max(trimode*10))/10,0.01),
#      cex.lab=1.4,xlab="Growth Rate (1/h)",main="")
r=trimode
popsim=popsimdat(r,N,time,iter)
popdata=popsim$PopData
poprate=popsim$SimRates
plotpopsimdat(popdata,r,time,iter)

par(op)