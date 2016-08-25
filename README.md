___

# MSc Submission: Modeling the consequences of heterogeneity in *S. cerevisiae* population dynamics.  

___

The current project state includes two available packages written along with the MSc Dissertation and all soure code required to generate the results presented in the submitted Dissertation. 

The python `muqfatc` package was written for analysing micro Quantiative Fitness Analysis images obtained by automated microscopy in order to obtain isogenic single lineage time course. 

The `detstochgrowth` package was written for subsequent time course analyses. Provided functions are for generating single lineage growth curves, generating population inferences from single lineage data and more. This package was used to analyze both the in-house muQFA data and the published data kindly provided to us by [Levy et al. 2012] (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3348152/). Synthetic data was also analyzed for exploratory purposes. These data are included in the Analyses/Synthetic directory. 

# Work Flow 

## Image Analysis

- Time course observations are saved as tiff images and are sorted into folders labeled by pin numbers (RxxCzz) specifying the rows and columns from the top left-hand corner of an agar plate on which yeast cells were grown.
- The image analysis tool `muqfatc` iterates through the pin directories analyzing single lineages observed in each time course. 
- Single blobs are captured by shape detection using a Canny edge map detection algorithm and optical dilation and erosion. Colonies are checked for circularity. 

*Insert Example Code Here*

- In order to be able to compare the population simulation done using single lineage data to true population simulations pin population growth can also be captured using `muqfatc`.

*Insert Example Code Here*

## Parameter Inference of Single Lineage Growth Curves 

Two frequentist, deterministic approaches are implemented. 

#### Log-linear regression

*Insert Example Code Here*

#### Exponential growth on the original scale 

*Insert Example Code Here*

A bayesian approach using `rjags` with various model specifications is available.

*Insert Example Code Here*

#### Exponential growth

*Insert Example Code Here*

#### Logistic growth 

*Insert Example Code Here*

#### Log-linear regression 

*Insert Example Code Here*

#### Why the B?
Both the exponential and the logistic model additional have a model specification with a B value. The B values sets a distinction between dividing and non-dividing growth curves such that N(t)=B*N0*exp(r*t) for the exponential model. B=1 means the lineage is growing / cells are dividing. B=0 means no growth. 
This was used to test as a possible option for the stochastic model outlined below. 

A Bayesian stochastic/deterministic hybrid model was also tested for parameter inference. NB: the current implementation is not suitable for a high-throughput level but can be used on individual growth curves which contain some noise and are relatively fast-growing. 

#### Discrete, stochastic birth-only growth 

*Insert Example Code Here*

## Population Simulations 

Population growth is simulated from single lineage data associated with a strain as follows:

*Insert Example Code Here*

Additional features which can be used to analyse the population simulations include the 

*Insert Example Code Here*

which assess lag duration with increasing inoculation size, growth rate with increasing inoculation sizes and single lineage percentile composition of population over time. 

## Comparing Population Simulations to Population Observations

*Insert Example Code Here*

## Image Analysis

- There is an additional feature in the `muqfatc` package which colour-codes the colonies on the plate according the growth rate and residual size. Code for this is included here but is not addressed in the Dissertation. No significant patters were found. 

____

# discstoch
Discrete stochastic kinetics for individual microbial microcolony growth curves.

![Yeast growing on agar surface](https://farm8.staticflickr.com/7169/6487645733_25284ce92d_z.jpg "Yeast microcolonies growing on agar")

Some [slides](http://lwlss.net/talks/uqfa) about microQFA.

Some [slides](http://lwlss.net/talks/discstoch) about discrete stochastic modelling of microcolony growth curves.

Data obtained from [Levy et al. 2012] (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3348152/) and [Ziv et al. 2013] (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3840306/).



