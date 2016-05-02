# Practicing a Linear Regression in RJAGS
# http://www.johnmyleswhite.com/notebook/2010/08/20/using-jags-in-r-with-the-rjags-package/

library(rjags)

N=1000
x=1:N
epsilon=rnorm(N, 0, 1)
y=x+epsilon 
#y is the data; numbers 1 to 1000 with some normally distributed added errors

write.table(data.frame(X = x, Y = y, Epsilon = epsilon),
            file = 'RJAGS.example.data',
            row.names = FALSE,
            col.names = TRUE)

View(read.table('RJAGS.example.data',header=TRUE))

#Bayesian model for regression written in RJAGS_Example.bug (Model specification file)
# model{                            ---->  have to tell JAGS that it's a model
#   for (i in 1:N){                 ---->  setting up the model for every single data point
#     y[i] ~ dnorm(y.hat[i],tau)    ---->  y[i] is distributed normally with mean y.hat[i] and precision tau
#                                          where precision is the reciprocal of the variance 
#     y.hat[i] <- a + b*x[i]
#   }
# a ~ dnorm(0, .0001)               ---->  non-informative normal priors assigned to a (mean 0, standard deviation 100)
# b ~ dnorm(0, .0001)               ---->  non-informative normal priors assigned to b (mean 0, standard deviation 100)
# tau <- pow(sigma, -2)
# sigma ~ dunif(0, 100)             ---->  non-informative uniform prior over the interval [0,100] assigned to the standard deviation sigma
#                                          which is deterministically transformed into tau
#                                          say that tau is a deterministic function (<- instead of ~) of sigma, after raising sigma to the power -2
# }

jags<-jags.model('RJAGS_Example.bug',
                data=list('x'=x, 'y'=y, 'N'=N), #names must be those in the JAGS model specification
                n.chains=4, #how many parallele chains to run
                n.adapt=100) #how many samples to throw away as part of the adaptive sampling period for each chain
update(jags,1000)
jags.samples(jags,c('a','b','tau'),1000) #draw 1000 sample from the sampler for the values of the named variables a and b 
#e.g. a=-0.1032664
#e.g. b=1.000242
#e.g. tau=1.005568