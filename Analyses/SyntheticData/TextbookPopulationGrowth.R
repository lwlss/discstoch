# This file generates a figure of what is considered typical population growth behaviour
# Used as Figure 1 in Dissertation

# Logistic growth model
logistic = function(K,r,x0,t) (K*x0*exp(r*t))/(K+x0*(exp(r*t)-1))

# Logistic growth model with lag (add delay)
lag_logistic = function(K,r,x0,lag,t) pmax(x0,logistic(K,r,x0,t-lag))

# Full growth curve (lag phase, exponential phase, stationary phase, decay phase)
full_curve = function(K,r,x0,lag,x0_d,d,t) lag_logistic(K,r,x0,lag,t)-lag_logistic(K,d,x0_d,lag,t)

op=par(mfrow=c(1,2))
# Original Scale
curve(full_curve(500,0.1,0.0001,30,0.00000001,0.05,x),
      from=0,to=700,log="",xlab="Time (h)",ylab="No. of Replicating Cells",lwd=4,cex.lab=1.4)

# Log Scale
curve(full_curve(500,0.1,0.0001,30,0.00000001,0.05,x),
      from=0,to=700,log="y",xlab="Time (h)",ylab="No. of Replicating Cells",lwd=4,cex.lab=1.4)
par(op)
