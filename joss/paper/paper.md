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
generally, for simulating from high-dimensional probability
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
language for specifying general hierarchical model structures, and
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
valid assignment of samplers to each unobserved model
dimension, mixing and matching between those samplers provided with
`nimble`. The base `nimble` package provides a variety of non-derivative-based samplers, including random walk Metropolis-Hastings [@robert1999metropolis],
slice sampling [@neal2003slice], conjugate samplers
[@george1993conjugate], and many others.  After configuration is finished, an MCMC
algorithm is generated according to the sampler assignments therein,
and executed to generate a sequence of samples.

The `nimbleHMC` package provides an implementation of HMC sampling which
is compatible for use within `nimble`.  Specifically, `nimbleHMC`
implements the No-U-Turn variety of HMC [HMC-NUTS; @hoffman2014no],
which removes the necessity of hand-specifying tuning parameters of
the HMC sampler.  Using
`nimbleHMC`, HMC samplers can be assigned to any subset of
continuous-valued model dimensions at the time of `nimble`'s MCMC
configuration, which may be used in combination with any other
samplers provided with the base `nimble` package.



# Statement of need

HMC is recognized as a state-of-the-art MCMC sampling strategy,
capable of efficiently generating samples with strong inferential
power.  A testimony to this, packages such as
`Stan` have built software exclusively around the use of HMC sampling.
As a result, however, such software is unable to operate on models with discrete (non-continuous)
valued dimensions, a result of the non-applicability of HMC.  Models with
discrete-valued dimensions arise in a range of common statistical motifs
including hidden Markov models, finite mixture models, and generally in
the presence of unobserved categorical data, among others
[@bartolucci2022discrete].  In contrast, other mainstream MCMC packages such as
`WinBUGS`, `OpenBUGS` and `jags` have the ability to sample discrete model dimensions,
but do not implement HMC.  This leaves a gap, as there is no support
for applying HMC sampling to continuous-valued dimensions of
hierarchical models which also contain discrete dimensions.

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
important from both practical and theoretical standpoints.
Indeed, there even exist packages such as `compareMCMCs`
[@de2022comparemcmcs], the purpose of which is to compare the relative
performance of distinct MCMC algorithms.

The `nimble` package uniquely provides the ability to custom-specify the
assignment of MCMC sampling algorithms, which allows the exploration and the
study of efficient approaches to MCMC.  No existing software to date
can operate on discrete model dimensions *and* offers the option of
HMC sampling.  Here, the `nimbleHMC` package augments the `nimble`
package by providing an HMC sampler suitable for use within `nimble`'s
MCMC.  This fills the gap, allowing HMC samplers to operate alongside the
existing continuous and discrete sampling algorithms available in `nimble`.





<!--
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

For a quick reference, the following citation commands can be used: 
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 
-et al., 2002)"



# Figures

 Figures can be included like this: 
![Caption for example figure.\label{fig:example}](figure.png) 
and referenced from text using \autoref{fig:example}. 

Figure sizes can be customized by adding an optional second parameter: 
![Caption for example figure.](figure.png){ width=20% }


# Acknowledgements
We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.
-->

# References

