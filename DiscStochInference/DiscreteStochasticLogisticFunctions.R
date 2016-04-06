# Discrete stochastic model simulating logistic growth (self-limiting population dynamics)
# Cell + Nutrients -> 2*Cell
# Mass-action kinetics, with rate constant r, carrying capacity K, initial population size c0

simCells=function(K,r,N0){
	# Unusually, for this model, we know the number of events a priori
	eventNo=K-N0
	# So we can just generate all required random numbers (quickly) in one go
	unifs=runif(eventNo)
	# Every event produces one cell and consumes one unit of nutrients
	clist=(N0+1):K
	nlist=K+1-clist
	# Simulate time between events
	dts=-log(1-unifs)/(r*clist*(1-clist/K))
	return(data.frame(t=c(0,cumsum(dts)),c=c(N0,clist)))
}

detLog=function(K,r,c0,t){
	return(K*c0*exp(r*t)/(K+c0*(exp(r*t)-1)))
}

simCellsHybrid=function(K,r,N0,NSwitch,detpts=100){
	# Every event produces one cell and consumes one unit of nutrients
	if(NSwitch>N0){
		# Unusually, for this model, we know the number of events a priori
		eventNo=NSwitch-N0
		# So we can just generate all required random numbers (quickly) in one go
		unifs=runif(eventNo)
		clist=(N0+1):NSwitch
		# Time between events
		dts=-log(1-unifs)/(r*clist*(1-clist/K))
		# Absolute times
		ats=cumsum(dts)
		tmax=max(ats)
	}else{
		clist=c()
		ats=c()
		tmax=0
	}
	
	# Switch to discrete deterministic logistic function
	clistdet=seq(NSwitch+(K-NSwitch)/detpts,K,(K-NSwitch)/detpts)
	tsdet=log((clistdet*(K - NSwitch))/((K - clistdet)*NSwitch))/r
	return(data.frame(t=c(0,ats,tmax+tsdet),c=c(N0,c(clist,clistdet))))
}

# Generate distriubtion of time to K-epsilon
getEnds=function(N0,Reps,K,r){
	ends=replicate(Reps,simCells(K,r,N0)$t[K/2])
	return(ends)
}

# Hybrid model expressed as number of cells at time t1, after starting at t0
simDt=function(K=1000,r=1,N0=1,NSwitch=100,t0=0,t1=1){
	if(NSwitch>N0){
		# Unusually, for this model, we know the number of events a priori
		eventNo=NSwitch-N0
		# So we can just generate all required random numbers (quickly) in one go
		unifs=runif(eventNo)
		clist=seq(N0,NSwitch)
		# Time between events
		dts=-log(1-unifs)/(r*clist[seq(2,(eventNo+1))]*(1-clist[seq(2,(eventNo+1))]/K))
		# Absolute times
		ats=c(t0,t0+cumsum(dts))
		tmax=max(ats)
		if(tmax>=t1){
			# Interpolate for estimate of c at t1
			af=approxfun(ats,clist,method="constant")
			return(af(t1))	
		}else{
			# Deterministic simulation from tmax to t1
			return(detLog(K,r,NSwitch,t1-tmax))
		}
	}else{
		return(detLog(K,r,N0,t1-t0))
	}
}