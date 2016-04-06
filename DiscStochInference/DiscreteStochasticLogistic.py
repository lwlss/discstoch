import numpy
import matplotlib.pyplot as plt
# Mass action kinetics with rate constant r
# Cell + Nutrition -> 2*Cell

def simPop(K,r,c0):
    # K: carrying capacity, r: growth rate, c0: inoculum density
    # Initial conditions
    c,n,t,i=c0,K,0,0
    # Storage
    res=numpy.zeros((K,2),dtype=numpy.float)
    res[i]=[0.0,c]

    # Simulation
    while c<K:
        t=t+numpy.random.exponential(1.0/(r*c*n))
        c+=1
        n-=1
        i+=1
        res[i]=[t,c]
    return(res[0:(i+1),:])

def simPop(K,r,c0):
    # K: carrying capacity, r: growth rate, c0: inoculum density
    # Initial conditions
    c,n,t=c0,K,0
    # Storage
    res=numpy.zeros((K,3),dtype=numpy.float)
    res[:,0]=numpy.log(1-numpy.random.random(K))
    res[0,1:]=[0.0,c]
    # Simulation
    res[1:,1:]=numpy.transpose(numpy.vstack((res[0:(K-1),1]+r*res[0:(K-1),0]*res[0:(K-1),2]*(K-res[0:(K-1),2]),res[0:(K-1),2]-1)))
    return(res[0:(K-c0),1:])

for j in xrange(0,5):
    res=simPop(K=1000,r=1,c0=10)
    plt.plot(res[:,0],res[:,1])
plt.xlabel("Time")
plt.ylabel("Number of cells")
plt.show()

