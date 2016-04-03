#Stochastic, birth-only models of cell population growth 

Supervisor: Dr Conor Lawless (conor.lawless@ncl.ac.uk - Institute for Cell and Molecular
Biosciences)

MSc Student: Helena Herrmann (h.a.herrmann2@newcastle.ac.uk - School of Computing Science)

######Full project proposal soon to follow. 

Deterministic modelling of population growth typically involves learning about an overall
intrinsic growth rate, incorporating the net effect of birth processes (e.g. cell division) and
death processes (e.g. apoptosis). On the other hand, a great advantage of stochastic
modelling of population dynamics is the ability to treat and learn about birth and death
processes separately [1]. As a result, much of the effort that has gone into the theory of
stochastic population modelling has concentrated on capturing population fluctuations
resulting from the conflict between birth and death processes. However, there are some
important biological model systems for which the death processes is negligible in vitro but for
which discrete stochastic models are still useful. Two examples are the growth of microbial
colonies and the growth of healthy human fibroblast populations.

When introduced to a new environment with abundant resources, microbial populations are
thought to undergo a preparatory lag phase before initiating cell division and exponential
growth [2]. However we have made microscopic observations of the first few divisions of
clonal Saccharomyces cerevisiae colonies (colonies inoculated from single cells) showing
cells that begin to divide almost immediately, but where some clones divide more quickly
than others. In these cases, it is plausible that population-level observations of a lag phase
are in fact the result of competition between clonal cell lineages whose rate of cell division is
stochastic. Much of the delay apparent at the population level could simply be the time taken
for the faster growing colonies to out-compete their slower-growing neighbours, or for
colonies which happened to undergo their first few cell divisions more quickly to out-compete
their neighbours. There is some existing work on stochastic interpretation of lag phase, but it
could be expanded considerably [3].

You will develop continuous, discrete and hybrid birth-only stochastic models of population
growth. You will develop Bayesian model parameter inference workflows for learning about
microbial growth rates from microscopy timecourses. Subject to availability of funding, you
will undertake some microscopic observations of clonal lineages of Saccharomyces
cerevisiae cells to capture the lag & exponential phases of growth in great detail. You will
develop image analysis tools for the generation of growth curves from microscopy data and
use captured growth curves and the new models to infer growth rates for these strains. You
will package and document the models and inference workflows to maximise accessibility for
other researchers. You will write an article describing the development of the models and
their validation against microscopy data.

####References 

[1] [Biological applications of the theory of birth-and-death processes] (http://dx.doi.org/10.1093/bib/bbk006)

[2] [Lag Phase Is a Distinct Growth Phase That Prepares Bacteria for Exponential Growth and
Involves Transient Metal Accumulation.] (http://dx.doi.org/10.1128/JB.06112-11)

[3] [Comparison of Stochastic and Deterministic Concepts of Bacterial Lag] (http://dx.doi.org/10.1006/jtbi.1998.0673)


#Aims and Objectives

The majority of experimental techniques are not of sufficiently high resolution to examine the behaviour of individual cells.  To get around this, (micro-)biologists design experiments to examine pure, isogenic, clonal cell cultures instead.  Clonal cultures are lineages derived from a single cell.  By examining purified cell populations in order to learn about the behaviour of single cells, we assume that all members of such populations behave in the same way.  However, there is increasing evidence that there is considerable heterogeneity between members of isogenic cell populations.

In this project, we will develop models to capture observable levels of  growth rate heterogeneity in isogenic microbial cultures.  We will use these models to explore how intracellular variability can affect our interpretation of population growth rate observations.  In particular we will examine evidence for lag phase in microbial growth curves.

Data analyses and model exploration will be based on three available data set: (i) uQFA data produced by Lawless (and possiblible repeat of the experiment) (ii) Levy et al. 2012 and (iii) Ziv et al. 2013. All data sets are based on Saccharomyces cerevisiae growth. 

##Aims

######1. 	If possible, repeat the uQFA experiments with the aim of obtaining higher resolution data (depending on resource and time availability)

This will involve undertaking microscopic obervations of clonal *S. cerevisiae* cells to capture the lag and exponential phase of growth in great detail. Image analysis tools will be required for generating growth curves from the data.

######2.	Find a model which best captures heterogeneity in microbial growth curves that allows for a distinction between ineherent lineage stochasticities and population heterogeneity. 

This is based on the notion that isogenic cell growth exhibits intrinsic stochasticities which may drown in the noise of extrinsic heterogeneities when considering population dynamics. Addressing the effect of isogenic variation on population dynamics may yield further mechanistic insight for analysing population growth curves. 

######3.	Explore the implications which a stochastic model may have on the interpretation of various growth phases with a particular emphasis on the lag phase.  

With the lag phase, thus far, being the least explored growth phase biologically and mathematically, a modelled distinction between isogenic single lineage variation and population variation is likely to lead to new insights and mechanistic interpreations. It is suspected that when taking cell growth and division time into account, the lag phase may actually be a mere artefact of a wide growth rate distribution within the population; μQFA video data obtained by Lawless provides little evidence for a lag phase at first sight. Also, microbiologists often sample cell populations during the exponential phase in order to obtain better results. It is thus worth analysing how apparent population heterogeneity affects sampling from different phases.

## Objectives

1.	Preparation and Data Exploration 
    1. Image analysis and experimentation; try and obtain uQFA data with higher resolution. 
    2. Ensure data accessibility by converting data sets into a general format which can be accessed using R and Python.
    3. Write scripts to extract growth curves from each of the data sets. Look at pulling out single growth curves and pulling out all growth curves for a single spot or genotype. 
    4. Use calibration curves to convert area vs time curves to cell count vs time curves.

######Write up: Stage 1 – Introduction and Background Reading (reuse project proposal). Generate plots to visualise the raw data for why the impact of hetergoneity is being researched.

2.	Model Development 
    1. Develop parameter inference workflows to learn about microbial growth rates. 
    2. Fit a deterministic model to ~20 isogenic growth curves for each data set. 
    3. Fit a stochastic birth-only model to ~20 isogenic growth curves for each data set. 
    4. Fit a hybrid model to ~20 isogenic growth curves for each data set to find a potential trade-off between speed and accuracy. Can a sensible cut-off for model switching be determined? 
    6. Can we improve the models by adding more biological information? Consider using the time required for cell growth and division as a lower bound for time sampling in the algorithm of the stochastic model. Potentially consider more complicated models which include some kind of death, or alternatively a switch to a non-dividing state, while maintaining monotonic increasing growth curves.  Is there a low rate of death which can still give increasing curves? *If so a new project title would be in order*
    7. Package and document the models. 

######Write up: Stage 2 - Model development. State model assumptions, implications, and validity. Analyse how these model fits differ and why. generate plots to demonstrate the accuracy of each of the developed models. State each of the required parameters and the biological significance of the parameters in the context of the model.

3.	Analyzing model Fit
    1. Which model seems to exhibit the closet fit to the data. 
    2. Compare model fit using Bayes factor. Bayes factor naturally penalizes for over-fitting. Meta-priors may prove useful for hierarchical, numerical integration of all three data sets. 

######Write up: Stage 3 – Model Exploration and Validity

4.	Exploring the Mechanistic Implications 
    1. Explore heterogeneity in the data; analyse the growth rate distribution. 
    2. How much of the apparent heterogeneity in growth rate is reduced when considering stochasticity?
    3. Does a lag phase *per se* even exist within the data or can this be explained by heterogeneity?
    4. See how apparent population heterogeneity might be affected by sampling from different phases.

######Write up: Stage 4 – Mechanistic implications on lag phase. Include figures that emphasize each of the concluded implications.

5.	Putting it all together 

######Write up: Stage 5 – Discussion + putting it all together; final Dissertation for submission on August 26th

####Rough TimeLine

- April 18th – April 29th (2 weeks): **Stage 1**
- May 2nd – June 10th (6 weeks): **Stage 2**
- June 13th – July 8th (4 weeks): **Stage 3**
- July 11th – July 15th (1 week): Catch-up/ Holiday
- July 18th – July 29th (2 weeks): **Stage 4**
- August 1st – August 26th (4 weeks): **Stage 5** 


