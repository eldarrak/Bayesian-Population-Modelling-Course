# Bayesian Population Modelling Course
This folder contains scripts and data for my part of the Advanced Population and Community Ecology masters Course at RUG.

**Recommended reading:**
Kéry, Marc, and Michael Schaub. Bayesian population analysis using WinBUGS: a hierarchical perspective. Academic Press, 2011. [_eBook available to RuG students_](http://web.b.ebscohost.com.proxy-ub.rug.nl/ehost/ebookviewer/ebook?sid=291ddf41-1ec4-4676-8a0e-97eed4b6b07d%40pdc-v-sessmgr02&ppid=pp_Cover&vid=0&format=EB)

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
