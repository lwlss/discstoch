source("DiscreteStochasticLogisticFunctions.R")


# Simulate populations with carrying capacity of 1000 cells
Ksim=1e6
init=1000
rfast=0.35
rslow=0.14
tmax=70
ncolonies=50
png("SlowFast.png",width=1000,height=1000,pointsize=24)
plot(NULL,xlim=c(0,tmax),ylim=c(0,Ksim),xlab="Time (h)",ylab="Cell no.",main=paste(ncolonies,"clonal colonies"),cex.lab=1.5)
date()
for(i in 1:ncolonies){
	fast=simCells(Ksim,rfast,init)
	slow=simCells(Ksim,rslow,init)
	if(length(fast)>2000){
		fast=fast[seq(1, length(fast), length(fast)/100)]
		slow=slow[seq(1, length(slow), length(slow)/100)]
	}
	points(fast$t,fast$c,type="s",col="grey",lwd=2)
	points(slow$t,slow$c,type="s",col="pink",lwd=2)
}
date()
curve(detLog(Ksim,rfast,init,x),from=0,to=tmax,col="black",lwd=3,add=TRUE)
curve(detLog(Ksim,rslow,init,x),from=0,to=tmax,col="red",lwd=3,add=TRUE)
dev.off()

# Real QFA cultures start with ~100 cells and end with ~1.6M cells
# giving very deterministic-looking curves
QFA=simCells(1600000,4,100)
plot(QFA$t,QFA$c,type="s",xlab="Time",ylab="Cell no.")

simCellsHybrid(100,1,1,10)

# Simulating a realistic number of cells takes 0.25 seconds
multi=replicate(10,as.matrix(simCells(10000,4,100)))
system.time({replicate(10,simCells(1600000,4,100))})["elapsed"]/10
system.time({replicate(10000,detLog(1600000,4,100,(1:1100/1100)*50))})["elapsed"]/10000
system.time({replicate(10000,simCellsHybrid(1000000,4,1,1000)$t)})["elapsed"]/10000

date()
for (i in 1:100) tt=simCellsHybrid(1000000,4,1,1000)$t
date()

N0s=c(1:100,seq(100,1000,50))
#rs=seq(0.01,0.5,0.01)
endsn0=sapply(N0s,getEnds,Reps=1000,K=1e6,r=0.35)
#endsr=sapply(rs,getEnds,N0=1,Reps=500,K=16000)
boxplot(endsn0)
cvsn0=apply(endsn0,2,sd)/apply(endsn0,2,mean)
#cvsr=apply(endsr,2,sd)

svg("COV.svg")
plot(N0s,cvsn0*100,xlab="N0",ylim=c(0,max(cvsn0*100)),ylab="Coefficient of variation (%)",type="l",lwd=2,cex.lab=1.5)
abline(v=100,col="red",lwd=3,lty=2)
legend("topright",legend="QFA inoculum density",lwd=2,col="red",lty=2)
#plot(rs,cvsr,xlab="r",ylab="Standard Deviation",type="l",lwd=2)
dev.off()

# Simulate populations with carrying capacity of 1000 cells
Ksim=1e6
init=1000
rfast=0.35
rslow=0.14
tmax=70
ncolonies=50
png("COV.png",width=1000,height=2000,pointsize=24)
op=par(mfrow=c(2,1))
plot(N0s,cvsn0*100,xlab="N0",ylim=c(0,max(cvsn0*100)),ylab="Coefficient of variation (%)",type="l",lwd=2,cex.lab=1.5)
abline(v=100,col="red",lwd=3,lty=2)
legend("topright",legend="QFA inoculum density",lwd=2,col="red",lty=2)
#plot(rs,cvsr,xlab="r",ylab="Standard Deviation",type="l",lwd=2)

plot(NULL,xlim=c(0,tmax),ylim=c(0,Ksim),xlab="Time (h)",ylab="Cell no.",main=paste(ncolonies,"clonal colonies"),cex.lab=1.5)
date()
for(i in 1:ncolonies){
	fast=simCells(Ksim,rfast,init)
	slow=simCells(Ksim,rslow,init)
	if(length(fast)>2000){
		fast=fast[seq(1, length(fast), length(fast)/100)]
		slow=slow[seq(1, length(slow), length(slow)/100)]
	}
	points(fast$t,fast$c,type="s",col="grey",lwd=2)
	points(slow$t,slow$c,type="s",col="pink",lwd=2)
}
date()
curve(detLog(Ksim,rfast,init,x),from=0,to=tmax,col="black",lwd=3,add=TRUE)
curve(detLog(Ksim,rslow,init,x),from=0,to=tmax,col="red",lwd=3,add=TRUE)
par(op)
dev.off()

sdss=apply(endss,2,sd)
sdsf=apply(endsf,2,sd)
plot(sdss,xlab="N0",ylab="Standard Deviation",type="l",lwd=2)
points(sdsf,col="red",type="l")

K=5000
NSwitch=1000
r=1
N0=1
dt=0.1

plot(NULL,xlim=c(0,20),ylim=c(0,K),xlab="time",ylab="Cell No.")
for (z in 1:20){
	clist=c(N0)
	tlist=c(0)

	for (i in 1:500){
		newt=tlist[i]+dt
		newc=simDt(K,r,clist[i],NSwitch,tlist[i],newt)
		clist=c(clist,newc)
		tlist=c(tlist,newt)
	}
	points(tlist,clist,type="s")
}

curve(detLog(K,r,N0,x),from=0,to=20,col="red",add=TRUE,lwd=3)


