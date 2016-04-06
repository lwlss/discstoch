source("DiscreteStochasticLogisticFunctions.R")

noiseSD=20000
truK=1e6
trur=0.345
switchN=1000
truth=simCellsHybrid(truK,trur,1,switchN,500)
maxt=max(truth$t[truth$t<Inf])
truthapprox=approxfun(truth$t,truth$c)
tlist=seq(20,maxt,(maxt-20)/10)
data=data.frame(c=truthapprox(tlist),t=tlist)
data$c=data$c+rnorm(dim(data)[1],0,noiseSD)
rownames(data)=data$t
data$t=NULL

plot(truth$t,truth$c,type="s",xlab="Time (h)",ylab="Cell number",cex.lab=1.5)
points(rownames(data),data$c,col="red",pch=16)

stepSim=function(x0=1, t0=0, deltat=1, th = c(100,3))  simDt(th[1],th[2],x0,switchN,t0,t0+deltat)

dataLik<-function(x,t,y,log=TRUE,...){
	ll=sum(dnorm(y,x,noiseSD,log=TRUE))
	if(log)
		return(ll)
	else
		return(exp(ll))
}

simx0=function(N,t0,...){
	# Initial condition for all populations is 1 cell
	return(matrix(1,nrow=N))
}

pfMLLik=function (n, simx0, t0, stepFun, dataLik, data) 
{
    times = c(t0, as.numeric(rownames(data)))
    deltas = diff(times)
    return(function(...) {
        xmat = simx0(n, t0, ...)
	  w=matrix(nrow=n)
        ll = 0
        for (i in 1:length(deltas)) {
		# Replace apply function with for loop to avoid vectorising vectorised function
		# if statements, seq function and a:b notation don't play nice with vectorisation...
		for(j in 1:n) {
			xmat[j,]=stepFun(x0=xmat[j,],t0=times[i],deltat=deltas[i],...)
			w[j]=dataLik(xmat[j,],t=times[i+1],y=data[i,],log=FALSE,...)
		}
            if (max(w) < 1e-20) {
                warning("Particle filter bombed")
                return(-1e+99)
            }
            ll = ll + log(mean(w))
            rows = sample(1:n, n, replace = TRUE, prob = w)
		# Typecast to matrix here, as otherwise 1D matrix gets converted to list...
            xmat = as.matrix(xmat[rows, ])
        }
        ll
    })
}

mLLik=pfMLLik(10,simx0,0,stepSim,dataLik,data)

#print(mLLik())
print(mLLik(th=c(truK,trur)))

# MCMC algorithm
date()
iters=1000
tune=0.01
thin=100
params=c(K=truK,r=trur)
pmin=c(K=0.5*truK,r=0.5*trur)
pmax=c(K=2*truK,r=2*trur)
p=length(params)
ll=-1e99
parmat=matrix(0,nrow=iters,ncol=p)
colnames(parmat)=names(params)
# Main MCMC loop
for(i in 1:iters){
	if (i%%(iters/10)==0) message(paste(i,date()))
	for(j in 1:thin){
		paramsprop=pmin/2
		while(sum((paramsprop<pmin)|(paramsprop>pmax))>0){
			# Reject particles outside of range
			paramsprop=params*exp(rnorm(p,0,tune))
		}
		llprop=mLLik(paramsprop)
		if (log(runif(1)) < llprop-ll){
			params=paramsprop
			ll=llprop
		}
	}
	parmat[i,]=params
}
date()

mcmcSummary(parmat)
svg("Posterior.svg",width=7, height=21,pointsize=24)
op=par(mfrow=c(3,1))
plot(truth$t,truth$c,type="s",xlab="Time (h)",ylab="Cell number",cex.lab=1.5,lwd=3)
points(rownames(data),data$c,col="red",pch=16)

plot(density(parmat[,1]),xlab="K (cells)",main="",lwd=3,cex.lab=1.5)
points(density(exp(runif(1000,log(pmin[1]),log(pmax[1])))),type="l",lty=2)
abline(v=truK,col="red",lwd=3)
legend("topleft",legend=c("True value","Prior","Posterior"),lwd=c(3,1,3),col=c("red","black","black"),lty=c(1,2,1))

plot(density(parmat[,2]),xlab="r (1/h)",main="",lwd=3,cex.lab=1.5)
points(density(exp(runif(1000,log(pmin[2]),log(pmax[2])))),type="l",lty=2)
abline(v=trur,col="red",lwd=3)
par(op)
dev.off()

plot(density(exp(runif(100000,log(0.1*3),log(10*3)))),type="l",lty=2,xlim=c(1.5,15))


