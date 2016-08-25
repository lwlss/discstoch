## R package `detstochgrowth` version 0.1

Link: https://github.com/lwlss/discstoch

Contact: h.herrmann2@gmail.com 

This R package provides tools for analyzing single lineage microbial growth curves and was written for analyzing single 
lineage data provided by [Levy et al. 2012](http://journals.plos.org/plosbiology/article?id=10.1371%2Fjournal.pbio.1001325) 
as well as muQFA data. 

#### General Usage

The package can be divided into 4 main features:

1. Frequentist deterministic exponential modelling (`FreqDetExpMod`)

Associated functions allow for single lineage growth curves in a dataset to be subdivided according to genotype,
colony and identifier and to be modelled using either a log-linear regression model or the exponential model on the original scale.

2. Frequentist deterministic piece-wise regression modelling (`FreqDetRRMod`)

Associated functions can be used for simulating population growth from individual lineages such that: N(t)=\sum_{i=1}^{n}exp(r_{i}t).
Population growth can then be analysed using the functions specified in `PopSimAnalyses` which capture simulated population growth.
Lineage compositions for population growth as well as the effect of inoculation size on population growth over time can be assessed.

3. Bayesian deterministic modelling (`BayesDetMod`)

Here one main function 'bayesdet` can be used to specify growth curves on which to do Bayesian parameter inference using `rjags`.
Various types of models including the logistic and exponential model with normal and log-normal error measurement can be specified for inference.

4. Bayesian stochastic modelling (`BayesStochMod`): 

Associated function allow for Bayesian inference using a discrete/stochastic hybrid model. 
The model implementation was adapted from the exisiting `smfsb` package to fit a birth-only stochastic growth model
for early time points of single lineage observations. 

#### Example Usage 

Example data is yet to be included in the package. However, a full workflow for all function in action can be found in the 
Analysis/LineageData , Analyses/PopulationData and Analyses/SyntheticData folders on the Git Hub repository. 
Sample data is provided in all folders. 

#### Dependencies

All package dependencies have been specified and should be automatically installed with the package. Functions in this
package are dependent on the `bcp`, `segmented`, and `fishplot` packages as specified. 
