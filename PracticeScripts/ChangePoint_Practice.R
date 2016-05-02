# Exploring packages for change point analyses 

N=36
x=1:N
epsilon=rnorm(N, 0, 2)
x=x+epsilon 
b=rep(c(0.5,-0.5),each=18)
y=b*x+rnorm(36)

library(bcp) # Bayesian analysis of change point problems
op_bcp=bcp(y,x)
#linear regression assumptions: observations (x,y) are partitioned into blocks and that linear models are appropriate within each block

plot(op_bcp, main="Linear Regression Change Point Example")

#getting the break point and the probability of the breakpoint
max_prob=max(op_bcp$posterior.prob,na.rm=TRUE)
breakpoint_location=which(op_bcp$posterior.prob==max_prob)
print(max_prob)
print(breakpoint_location)

library(changepoint) #finds changes in the mean of the data 
#but I want to find changes in the slope of the data.... 

op_cp=cpt.mean(y,method="AMOC") #at most one change point
print(summary(op_cp))
plot(op_cp)

#definitely think the bcp provides a nicer solution than the changepoint package
#use that one to set a cutoff for when two use piece-wise regression and when to use to simpler model
#say if the max_prob for the breakpoint is greater than 0.5 ... 
#Except... changepoint package is much faster when dealing with say 1000 data points but both are fine for 30 data points