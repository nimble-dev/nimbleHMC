---
title: "nimbleHMC: An R package for Hamiltonian Monte Carlo sampling in nimble"
tags:
  - R
  - MCMC
  - HMC
  - nimble
authors:
  - name: Daniel Turek
    orcid: 0000-0002-1453-1908
    corresponding: true
    equal-contrib: true
    affiliation: 1
  - name: Perry de Valpine
    orcid: 0000-0002-8329-6796
    equal-contrib: true 
    affiliation: 2 
  - name: Christopher J. Paciorek
    equal-contrib: true 
    affiliation: 2 
affiliations:
 - name: Williams College, USA
   index: 1
 - name: University of California, USA
   index: 2
date: 13 August 2022
bibliography: paper.bib
output: pdf_document
---


<!--
setwd('~/github/nimble/nimbleHMC/joss/paper')
f <- 'paper.md'
library(rmarkdown)
rmarkdown::render(f, output_format = 'pdf_document')
system('open paper.pdf')
-->


# Summary

Markov chain Monte Carlo (MCMC) algorithms are widely used for
fitting hierarchical (graphical) models to observed data, and more
generally for simulating from high-dimensional probability
distributions.  MCMC is the predominant tool used in Bayesian
analyses, where the distribution of interest (the "target
distribution") is defined as the posterior
distribution of the unknown model parameters conditional on the data.
MCMC does not specify a single algorithm, but rather
a family of algorithms admitting any assignment of valid sampling
techniques ("samplers") to the unobserved model dimensions.
There exist a vast and diverse landscape of valid samplers to draw upon,
which differ significantly in their underlying approaches
to the sampling problem, complexity, autocorrelation of the samples
produced, and applicability.

Hamiltonian Monte Carlo [HMC; @brooks2011handbook] is one such sampling technique which can be
applied to any subset of continuous-valued model dimensions.  HMC uses
the gradient of the target distribution to generate large transitions
(in parameter space) in the output sequence of samples.  This results in low
autocorrelation, and therefore
high information content.  That is, the samples generated using HMC
are more likely to be highly informative about the target
distribution of interest, relative for example to an equal-length sequence of highly
autocorrelated samples.   This rich information content does not come
freely, however, as calculating gradients of the target distribution
is computationally expensive.

There exist numerous software packages
which provide implementations of MCMC for mainstream use such as
`nimble` [@de2017programming], `jags` [@plummer2003jags], `pyMC`
[@fonnesbeck2015pymc], and `Stan` [@carpenter2017stan].  Each such package provides a
language for specifying general higherarchical model structures, and
supplying data.  Following specification of the problem, each package generates an
MCMC algorithm which specifically samples from the target posterior
distribution of the specified model, and executes this algorithm to
generate a sequence of samples from this distribution.
These packages differ, however, in their approaches to sampler assignment for each
unobserved model dimension.  As sampling techniques vary in terms of
computational demands and the quality of the samples produced, the
effectiveness of the MCMC algorithms may vary depending on the
software used, and the particular model at hand.  Each software
package provides a valid, but distinct approach for assigning samplers
to define the MCMC algorithm.

Among general-purpose MCMC software packages, `nimble`
uniquely provides the ability to specify which samplers
are applied to each model dimension.  Prior to generating an
executable MCMC algorithm, `nimble` has the intermediate stage of MCMC
configuration.  At configuration time, users may select any
valid assignment of samplers to each unosberved model
dimension, mixing and matching between those samplers provided with
`nimble`. The base `nimble` package provides a variety of non-derivative-based samplers, including random walk Metropolis-Hastings [@robert1999metropolis],
slice sampling [@neal2003slice], conjugate samplers
[@george1993conjugate], and many others.  After configruation is finished, an MCMC
algorithm is generated according to the sampler assignments therein,
and executed to generate a sequence of samples..

The `nimbleHMC` package provides an implementation of HMC sampling which
is compatible for use within `nimble`.  Specifically, `nimbleHMC`
implements the No-U-Turn variety of HMC [HMC-NUTS; @hoffman2014no],
which removes the necessity of hand-specifying tuning parameters of
the HMC sampler.  Using
`nimbleHMC`, HMC samplers can be be assigned to any subset of
continuous-valued model dimensions at the time of `nimble`'s MCMC
configuration, which may be used in combination with any other
samplers provided with the base `nimble` package.




<!-- 
Markov chain Monte Carlo (MCMC) algorithms are used to simulate from  
complicated probability distributions. MCMC is very widely used to  
implement Bayesian statistical analysis, where the distribution of  
interest (the “target distribution”) is a posterior distribution of  
pa-rameters given data. In this context, the posterior is known only  
up to a constant, so only relative probabilities (or densities) can be  
easily calculated, which is sufficient for MCMC to work. In Bayesian  
statistical analysis, MCMC algorithms for large or complex statistical  
models and data sets are sometimes run for minutes, hours or days,  
making them an analysis bottleneck, so there is a premium on  
efficiency. MCMC efficiency includes both computational speed and  
algorithmic mixing, which refers to how well the algorithm explores  
the posterior distribution from one iteration to the  
next. Computational speed may comprise one or more steps such as  
algorithm setup, MCMC “burn-in” or “warm-up” phases, and MCMC  
execution or “sampling.”
 
There are many MCMC algorithms (also called “samplers”) and software  
packages implement-ing them. Because MCMC samplers can be validly  
combined (e.g., iterated in sequence), for example with different  
samplers for different dimensions of a target distribution, there is  
an enormous space of MCMC methods. Invention of new methods,  
comparisons among methods, and theoretical study of MCMC mixing are  
all important areas of active research. Various soft-ware packages  
provide samplers such as Gibbs, adaptive random-walk  
Metropolis-Hastings, slice, Hamiltonian, multivariate (“block”) or  
other variants of these, and others. Different MCMC algorithms can  
yield efficiencies that differ by orders of magnitude for a particular  
problem, with these variations in efficiency being problem-dependent.  
 
The R package compareMCMCs provides a highly modular system for  
managing performance comparisons among MCMC software packages for  
purposes of research on MCMC methods. MCMC runs can take a long time,  
so the output (“samples”) and components of compu-tation time from a  
run are stored regardless of whether performance metrics are computed  
immediately. Arbitrary MCMC packages (“MCMC engines”) can be added to  
the system by writing a simple plug-in or wrapper to manage inputs and  
outputs in a unified way. Con-versions among model parameter names  
and/or different parameterizations can be provided to standardize  
across packages. Performance metrics are organized by model parameter  
(one result per parameter per MCMC engine), by MCMC (one result per  
MCMC engine), or arbi-trarily (a user-defined list of metric results  
per MCMC engine). Built-in metrics include two methods of estimating  
effective sample size (ESS), posterior summaries such as mean and  
common quantiles, efficiency defined as ESS per computation time, rate  
defined as compu-tation time per ESS, and minimum efficiency per  
MCMC. New metrics can be provided by a plug-in system and applied  
programmatically to a set of MCMC samples without re-running the MCMC  
engines. Finally, standardized graphical comparison pages can be  
generated in html. Built-in graphical outputs include figures  
comparing MCMC efficiency and/or rate on a per-parameter or per-MCMC  
basis as well as comparing posterior distributions. New graphical  
outputs can be provided by a plug-in system. In summary, compareMCMCs  
is modular and extensible for running new MCMC engines on comparable  
problems, for creating new metrics of interest (e.g., posterior  
summaries or effective sample size estimated in different ways), and  
for creating new graphical comparison outputs into a report.  
 
Use of compareMCMCs supports but does not require a primary role for  
MCMCs created with the nimble (de Valpine et al., 2017, 2021) package  
for hierarchical statistical models. That is because nimble provides  
greater flexibility than other packages to customize its MCMC system,  
configuring which samplers will operate on which parts of a model  
and/or writing new samplers. Thus, it is of interest to compare  
multiple MCMC methods all implemented within nimble. Furthermore,  
nimble uses a model language that is a dialect of that used by
WinBUGS, OpenBUGS, MultiBUGS, and JAGS (Goudie et al., 2020; D. Lunn
et al., 2009; D. J. Lunn et al., 2000; Plummer & others, 2003). These
packages are often called from R via packages such as R2WinBUGS
(Sturtz et al., 2005), rjags (Plummer, 2019), and jagsUI (Kellner,
2019). Therefore, for fully compatible models, comparisons between
nimble and JAGS can be run in compareMCMCs from the same model and
data specifications. A plug-in is also provided for Stan via rstan
(Stan Development Team, 2020), and the extension system to plug in new
MCMC engines is clearly documented.
-->




# Statement of need

HMC is recognized as a state-of-the-art MCMC sampling strategy,
capable of efficiently generating samples with strong  inferential
power.  A testimony to this, packages such as
`Stan` have built software exclusively around the use of HMC sampling.
As a result, however, such software is unable to operate on models with discrete (non-continuous)
valued dimensions, a results of the non-applicability of HMC.  Models with
discrete-valued dimensions arise in a range of common statistical motifs
including hidden Markov models, finite mixture models, and generally in
the presence of unobserved categorical data, among others
[@bartolucci2022discrete].  In constrast, other mainstream MCMC packages such as
`WinBUGS`, `OpenBUGS` and `jags` have the ability to sample discrete model dimensions,
but do not implement HMC.  This leaves a gap, as there is no support
for applying HMC sampling to continuous-valued dimensions of
hierchical models which also contain discrete dimensions.

It is an open question what MCMC algorithm, or which combination of
samplers, will optimize the fitting of any particular hierarchical model and dataset.
The metric of interest is the effective sample size of the sequence of
samples (which accounts for autocorrelation) generated per unit
runtime of the MCMC.  That is, how quickly an MCMC algorithm generates
meaningful information to characterize the target distribution.  This
metric is called MCMC efficiency [@turek2017automated], but what assignment of samplers
maximizes this metric is a difficult and open question
[@ponisio2020one].  For that reason, the ability to mix-and-match
samplers from among as large a pool of candidates as possible is
important from both practical and a theoretical standpoints.
Indeed, there even exist packages such as `compareMCMCs`
[@de2022comparemcmcs], the purpose of which is to compare the relative
performance of distinct MCMC algorithms.

The `nimble` package uniquely provides the ability to custom-specify the
assignment of MCMC sampling algorithms, which allows the exploration and the
study of efficiency approachs to MCMC.  No existing software to date
can operate on discrete model dimensions *and*  offers the option of
HMC sampling.  Here, the `nimbleHMC` package augments the `nimble`
package by providing an HMC sampler suitable for use within `nimble`'s
MCMC.  This fills the gap, allowing HMC samplers to operate alongside the
existing continuous and discrete sampling algorithms avaialble in `nimble`.







<!--
Many other packages run MCMC algorithms and/or post-process MCMC
results, but compa reMCMCs is distinct in its goal of supporting MCMC
research by comparing MCMC methods. Packages that run MCMC from R are
documented on the “Cran Task View” page for “Bayesian Inference”
(Park, 2021) of the Comprehensive R Archive Network (CRAN). Some
popular general packages include those listed above as well as others
such as MCMCpack (Martin et al., 2011) and LaplacesDemon (Statisticat
& LLC., 2021). Furthermore, there are MCMC engines based in Python,
such as PyMC (Salvatier et al., 2016), and other languages. These may
be called via appropriate interfaces from R to other languages.

Of the packages listed on the “Bayesian Inference” Task View, only the
SamplerCompare (Thompson, 2011) package appears to specifically
support the goal of comparing MCMC performance. However, this package
can only compare MCMC samplers that have exactly one scalar tuning
parameter, target distributions that are continuous with constant
dimension, and are implemented within the package.

Packages for post-processing of MCMC samples (e.g., coda (Plummer et
al., 2006), BayesP ostEst (Scogin et al., 2019), and MCMCvis
(Youngflesh, 2018)) aim to provide features for scientific summary and
presentation of results, whereas compareMCMCs provides features for
comparisons of algorithm performance across packages. Assessing MCMC
performance is not simply a matter of computational benchmarking. For
example, effective sample size is itself a non-trivial property to
estimate by statistical methods, different metrics may be of interest
for different purposes, and consistency of algorithm results between
different MCMC engines can only be determined statistically,
i.e. within simulation error. Therefore, the features needed for
comparing MCMC performance are distinct from those needed for
presenting scientific results based on MCMC.
-->


# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

<!--
For a quick reference, the following citation commands can be used: 
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2
-et al., 2002)"
-->


# Figures

<!--
Figures can be included like this: 
![Caption for example figure.\label{fig:example}](figure.png) 
and referenced from text using \autoref{fig:example}. 

Figure sizes can be customized by adding an optional second parameter: 
![Caption for example figure.](figure.png){ width=20% }
-->


# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References

