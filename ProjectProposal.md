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
