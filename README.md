# Bayesian Population Modelling Course
This folder contains scripts and data for my part of the Advanced Population and Community Ecology masters course at the University og Groningen.

**Recommended reading:**
Kéry, Marc, and Michael Schaub. Bayesian population analysis using WinBUGS: a hierarchical perspective. Academic Press, 2011. [_eBook available to RuG students_](http://search.ebscohost.com.proxy-ub.rug.nl/login.aspx?direct=true&db=nlebk&AN=407875&site=ehost-live&scope=site&ebv=EB&ppid=pp_iii)

The course is based on Just another Gibbs Sampler (or shortly JAGS) software.
JAGS code for all exercises in the book is available here https://www.vogelwarte.ch/de/projekte/publikationen/bpa/code-for-running-bpa-using-jags

Before start of the course, **please install all the required software**:

1. install R (you may also want to RStudio, notepad++ or other GUI of your choise);
2. install JAGS http://mcmc-jags.sourceforge.net;
3. install additional R packages that we will use during the course:
```{r}
install.packages ('lme4')
install.packages ('lattice')
install.packages ('coda')
install.packages ('R2WinBUGS')
install.packages('rjags')
install.packages ('R2jags')
```
You can also source some useful functions we will be useing throughout the course
```{r}
source('https://git.io/fpufJ')
```

1. GLM in R and JAGS [link to the script](https://github.com/eldarrak/Bayesian-Population-Modelling-Course/blob/master/1-GLMS.R)
2. Mixed effects models in JAGS [link to the script](https://github.com/eldarrak/Bayesian-Population-Modelling-Course/blob/master/2-Mixed_effects_GLMs.R)
3. State-space models [link to the script](https://github.com/eldarrak/Bayesian-Population-Modelling-Course/blob/master/3-state-space-models.R)
4. Capture-recapture models 1 [link to the script](https://github.com/eldarrak/Bayesian-Population-Modelling-Course/blob/master/4-capture-recapture-1.R)
5. Real data analyses [link to the script](https://github.com/eldarrak/Bayesian-Population-Modelling-Course/blob/master/real-data_analyses.R)

Note that current course havily relies on the Bayesian population analysis book but does not cover it completely, so please read the book if you want to go further.

Capture-recapture modelling can also be done in the frequentist framwork. The main software for this analysis is [program Mark](http://www.phidot.org/software/mark/) by Gary White. [RMark](https://cran.r-project.org/web/packages/RMark/index.html) R package is an interface to Mark from R written by Jeff Laake. The great explanation of all conemporary CR methods by Evan Cooch is available [here](http://www.phidot.org/software/mark/docs/book/). 

If you run in to problems with capture-recapture methods ask you question [here](http://www.phidot.org/forum/index.php).


