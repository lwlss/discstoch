___

# Modeling the consequences of heterogeneity in *S. cerevisiae* population dynamics.  

___

The current project state includes two available packages written along with the MSc Dissertation and all soure code required to generate the results presented in the submitted Dissertation. 

The python `muqfatc` package was written for analysing micro Quantiative Fitness Analysis images obtained by automated microscopy in order to obtain isogenic single lineage time course. 

The `detstochgrowth` package was written for subsequent time course analyses. Provided functions are for generating single lineage growth curves, generating population inferences from single lineage data and more. This package was used to analyze both the in-house muQFA data and the published data kindly provided to us by [Levy et al. 2012] (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3348152/). Synthetic data was also analyzed for exploratory purposes. These data are included in the Analyses/Synthetic directory. 

# Work Flow 

## 1) Image Analysis

- Time course observations are saved as tiff images and are sorted into folders labeled by pin numbers (RxxCzz) specifying the rows and columns from the top left-hand corner of an agar plate on which yeast cells were grown.
- The image analysis tool `muqfatc` iterates through the pin directories analyzing single lineages observed in each time course. 
- Single blobs are captured by shape detection using a Canny edge map detection algorithm and optical dilation and erosion. Colonies are checked for circularity. 

```
FINALPHOT=20 # Time course length
fullpath = os.path.abspath(syspath) # Working directory

# Generating a grey background of median pixels of all first images
DX,DY =20,20 # Background size
bk=muqfatc.makeBackground(myfolders,DX,DY)   

# Generating time course estimates (as output matrices saved in text files)
# and time course images (saved in respective output folders). 
muqfatc.lineagetimecourse(myfolders,FINALPHOT,bk,DX,DY,fullpath,apply_filt=True,
                       save_pics=True,log=False,show_im=False)
```

- In order to be able to compare the population simulation done using single lineage data to true population simulations pin population growth can also be captured using `muqfatc`.

```
# Getting time and area for pin growth estimates for 90 observations
area,time=pinia.pintimecourse(myfolders,90,fullpath)
```

## 2) Parameter Inference of Single Lineage Growth Curves 

Two frequentist, deterministic approaches are implemented. 

#### Log-linear regression (Frequentist)

```
rates=c()
for (i in 1:dim(area)[1]){
  k=detstocgrowth::LM_growthrate(area[i,],times[i,])$rate
  rates=c(rates,k)
}
```

#### Exponential growth on the original scale (Frequentist)

```
rates=c()
for (i in 1:dim(area)[1]){
  k=detstocgrowth::EXP_growthrate(area[i,],times[i])$rate
  rates=c(rates,k)
}
```

A bayesian approach using `rjags` with various model specifications is available. A few examples are shown below.

#### Exponential growth (Bayesian)

```
gc=250:255 # Inference for growth curves 250:255
area=area[gc,]
times=times[gc,]
data=data[gc,]

iterations=1000000
thin=1000
no_chains=4
model_choice=2
bayesinf=detstocgrowth::BayesDet(area,times,model_choice,iterations,thin,no_chains)
```

#### Logistic growth (Bayesian)

```
model_choice=1
bayesinf=detstochgrowth::BayesDet(area,times,model_choice,iterations,thin,no_chains)
```

#### Log-linear regression (Bayesian)

```
model_choice=3
bayesinf=detstochgrowth::BayesDet(area,times,model_choice,iterations,thin,no_chains)
```

#### Why the B?
Both the exponential and the logistic model additional have a model specification with a B value. The B values sets a distinction between dividing and non-dividing growth curves such that N(t)=B*N0*exp(r*t) for the exponential model. B=1 means the lineage is growing / cells are dividing. B=0 means no growth. 
This was used to test as a possible option for the stochastic model outlined below. 

A Bayesian stochastic/deterministic hybrid model was also tested for parameter inference. NB: the current implementation is not suitable for a high-throughput level but can be used on individual growth curves which contain some noise and are relatively fast-growing. 

#### Discrete, stochastic birth-only growth (Bayesian)

```
# Apply calibration and format data 
cells=as.numeric(area)/area[,1]
data=data.frame(c=cells,t=as.numeric(times))
rownames(data)=data$t
data$t=NULL

# Specify parameters for Inference
guessK=100000
guessr=0.3
params=c(K=guessK,r=guessr)
pmin=c(K=0.5*guessK,r=0.01)
pmax=c(K=2*guessK,r=1)
switchN=100 # switch from stochastic to deterministic model
noiseSD=20
iterations=1000000
thin=10000
tune=0.02

# Inference using the Hybrid Model
detstochgrowth::mcmc_chain=BayesHybrid(params,pmin,pmax,iterations,tune,thin,noiseSD,switchN,data)
smfsb::mcmcSummary(mcmc_chain,show=TRUE,plot=TRUE)

# Plotting the MCMC Output
detstochgrowth::PlotPosteriorProb(pmin,pmax,mcmc_chain)

# Plotting the Posterior Predictives
detstocgrowth::PlotPosteriorPredictive(data,mcmc_chain,switchN,guessK)
```

## 3) Population Simulations 

Population growth is simulated from single lineage data associated with a strain as follows:

```
strain_names=unique(data$genotype)
pickstrain="HIS3"
strain=detstocgrowth::subset_strain(data,area,times,pickstrain)
strain_rates=rates[which(data$genotype==pickstrain)]

iterations=100
yl=10^22
gr=seq(0,0.5,0.005)
N=1000
t=seq(1,48)
gr=seq(0,0.5,0.001)
simpopdat=detstocgrowth::pop_sim_plot(strain,strain_rates,gr,yl,iterations,t)
```

Additional features which can be used to analyse the population simulations include

```
# Generating a colour map according to which population simulations are coloured
detstochgrowth::yellowredmap(gr)

# Producing fishplots
N=10000
sample_rates=sample(strain_rates,N,replace=TRUE)
time=seq(1,48)
detstochgrowth::pop_comoposition(sample_rates,time,gr)

# Calculating percentile lineage contributions
t=seq(1,48,1)
iter=1000
StartPops=c(50,100,500,1000,5000,10000)
#StartPops=c(50,100,500)
for (N in StartPops){
  newinfo=detstochgrowth::pop_composition_dat(t,iter,N,strain_rates)
  if (N==StartPops[1]){
    info=newinfo
  } else{
    info=data.frame(info,newinfo)
  }
}

# Estimating lag phase durations
N=seq(100,10000,100)
time=seq(0,48,1)
iterations=100
strainlag=detstocgrowth::lagduration(strain=0,strain_rates,t,N,iterations)
```

which assess lag duration with increasing inoculation size, growth rate with increasing inoculation sizes and single lineage percentile composition of population over time. 

## 4) Comparing Population Simulations to Population Observations

```
# Plotting the Pin Population Estimates
PinDat=detstochgrowth::pin_pop_plot(strain,strainarea,sraintimes)

sortedfolders=sortedfolders[18] #Looking only at pin 18
obs=30 # Time course length
conv=92.5 # Conversion factor; area/conv=cells
it=100 #Number of simulations
Comp=detstocgrowth::PlotPinObsSim(area,times,sortedfolders,obs,conv,it,leg=FALSE)
```

## 5) Image Analysis

- There is an additional feature in the `muqfatc` package which colour-codes the colonies on the plate according the growth rate and residual size. Code for this is included in the Analyses/ImageAnalysis directory but is not addressed in the Dissertation. No significant patters were found. 

____

# discstoch
Discrete stochastic kinetics for individual microbial microcolony growth curves.

![Yeast growing on agar surface](https://farm8.staticflickr.com/7169/6487645733_25284ce92d_z.jpg "Yeast microcolonies growing on agar")

Some [slides](http://lwlss.net/talks/uqfa) about microQFA.

Some [slides](http://lwlss.net/talks/discstoch) about discrete stochastic modelling of microcolony growth curves.

Data obtained from [Levy et al. 2012] (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3348152/) and [Ziv et al. 2013] (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3840306/).



