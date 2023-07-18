---
title: "nimbleHMC: An R package for Hamiltonian Monte Carlo sampling in nimble"
tags:
  - R
  - hierarchical model
  - Markov chain Monte Carlo
  - Hamiltonian Monte Carlo
  - nimble
authors:
  - name: Daniel Turek
    orcid: 0000-0002-1453-1908
    corresponding: true
    equal-contrib: false
    affiliation: 1
  - name: Perry de Valpine
    orcid: 0000-0002-8329-6796
    equal-contrib: false
    affiliation: 2 
  - name: Christopher J. Paciorek
    equal-contrib: false
    affiliation: 2 
affiliations:
 - name: Lafayette College, USA
   index: 1
 - name: University of California, USA
   index: 2
date: 18 July 2023
bibliography: paper.bib
output: pdf_document
---






# Summary

Markov chain Monte Carlo (MCMC) algorithms are widely used for
fitting hierarchical models to data.
MCMC is the predominant tool used in Bayesian
analyses to generate samples from the posterior
distribution of model parameters conditional on observed data.
MCMC is not a single algorithm, but actually a framework which admits any assignment of sampling
techniques to unobserved parameters.
There exists a vast set of valid samplers to draw upon,
which differ in complexity, autocorrelation of the samples
produced, and applicability.  Hamiltonian Monte Carlo [HMC; @brooks2011handbook] is one such
sampling technique applicable to any subset of continuous-valued parameters.  HMC uses
the gradients to generate large transitions in the output sequence of samples.  This results in low
autocorrelation, and therefore the samples generated using HMC
are more likely to be highly informative about the target
distribution, relative for example to an equal-length sequence of highly
autocorrelated samples.   This rich information content does not come
freely, as calculating gradients is computationally expensive.

Many software packages offer implementations of MCMC, such as
`nimble` [@de2017programming], `WinBUGS` [@lunn2000winbugs], `jags` [@plummer2003jags], `pyMC`
[@fonnesbeck2015pymc], and `Stan` [@carpenter2017stan].
These packages differ, however, in their approaches to sampler
assignment to model parameters.  As sampling techniques vary in
computational demands and quality of the samples, the
effectiveness of the MCMC algorithms will vary depending on the
software and model; each software
package provides distinct approach for sampler assignment.

Among MCMC software packages, `nimble` uniquely allows specification
of which samplers are used. Users may select any
valid assignment of samplers to each parameter,
selecting from among the suite of non-derivative-based samplers
provided with `nimble`.  These include
random walk Metropolis-Hastings sampling [@robert1999metropolis],
slice sampling [@neal2003slice], elliptical slice sampling
[@murray2010elliptical], automated factor slice sampling [@tibbits2014automated],
conjugate sampling [@george1993conjugate], and many others.  The `nimbleHMC` package provides an implementation of HMC sampling
for use within `nimble`.  Specifically, `nimbleHMC`
implements the No-U-Turn variety of HMC [HMC-NUTS; @hoffman2014no].
HMC samplers can be assigned to any set of
continuous-valued parameters, and may be used in combination with other
samplers provided with `nimble`.





# Example

The following example demonstrates fitting a hierarchical
model to data using `nimbleHMC`.
We use the European Dipper \emph{Cinclus cinclus)} dataset drawn from ecological
capture-recapture
[\emph{e.g.}, @lebreton1992modeling; @turek2016efficient].
Modelling includes both continuous parameters to undergo HMC
sampling and discrete parameters which cannot be sampled via HMC.
This combination is not supported by software other than `nimbleHMC`.

Individual birds are captured, tagged, and potentially recaptured on
subsequent sighting occasions.  Data is a $294 \times 7$
binary-valued array of capture histories of 294 uniquely tagged birds
over 7 years.  Model parameters are detection
probability $p$, and annual survival rates on non-flood years $\phi_1$
and flood years $\phi_2$.  Data is provided in the R package `mra` [@mcdonald2018mra].

```
library(mra) 
data(dipper.data) 
y <- dipper.data[,1:7]
```

We specify the hierarchical model using uniform priors on the interval
$[0,1]$ for all parameters.
Binary-valued latent states $x_{i,t}$ represent the true alive (1) or dead (0)
state of individual $i$ on year $t$.  Doing so allows the survival
process to be modelled as $x_{i,t+1}~\sim~\text{Bernoulli}(\phi_{f_t} \cdot
x_{i,t})$ where $f_t$ indicates the flood/non-flood history of year
$t$, and observations are modelled as $y_{i,t}~\sim~\text{Bernoulli}(p \cdot x_{i,t})$.

```
library(nimbleHMC) 

code <- nimbleCode({
    phi[1] ~ dunif(0, 1)
    phi[2] ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for(i in 1:N) {
        for(t in (first[i]+1):T) {
            x[i,t] ~ dbern(phi[f[t]] * x[i,t-1])
            y[i,t] ~ dbern(p * x[i,t])
        }
    }
})
```

A `nimble` model object is now built.
The argument `buildDerivs = TRUE` affects derivatives of likelihood calculations
to be built into the model object to support
derivative-based algorithms -- here, HMC sampling.

```
Rmodel <- nimbleModel(
    code,
    constants = list(N = nrow(y), T = ncol(y), f = c(1,2,2,1,1,1,1),
                     first = apply(y, 1, which.max)),
    data = list(y = y),
    inits = list(phi = c(0.5, 0.5), p = 0.5, x = array(1, dim(y))),
    buildDerivs = TRUE)
```

Next we create an MCMC configuration object, which specifies the
sampling algorithm to be applied to each parameter.
By default, `configureMCMC` uses `nimble`'s default sampler
assignments of adaptive random walk Metropolis-Hastings
[`RW` sampler; @robert1999metropolis] for each parameter, and a
`binary` Gibbs sampler for each $x_{i,t}$ latent state.

```
conf <- configureMCMC(Rmodel)

## RW sampler (3)
##   - phi[]  (2 elements)
##   - p
## binary sampler (848)
##   - x[]  (848 elements)
```

Now we customize the MCMC configuration object to use HMC sampling for
the model parameters.  
`replaceSamplers` replaces current samplers operating on
$\phi_1$, $\phi_2$ and $p$ instead with the `HMC` sampler provided in
 `nimbleHMC`.

```
conf$replaceSamplers(target = c("phi", "p"), type = "HMC")
conf$printSamplers(byType = TRUE)

## HMC sampler (1)
##   - phi, p 
## binary sampler (848)
##   - x[]  (848 elements)
```

Alternatively, the convenience function `configureHMC(Rmodel)`
may be used to create an identical MCMC configuration, applying
HMC sampling to $\phi_1$, $\phi_2$ and $p$, and default binary
samplers for discrete parameters.

Now we build and compile the MCMC algorithm.

```
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
```

We execute the MCMC for 20,000 iterations, and discard the initial 10,000 samples
as burn-in.

```
set.seed(0)
samples <- runMCMC(Cmcmc, niter = 20000, nburnin = 10000)

##  [Note] HMC sampler (nodes: phi[1], phi[2], p) is using 1000 warmup iterations. 
##  [Note] HMC sampler (nodes: phi[1], phi[2], p) encountered 3 divergent paths. 
```

The HMC sampler outputs two notes, indicating the number
of warmup iterations and the total number of divergent paths
encountered [see @hoffman2014no for details].

Finally, posterior summary statistics are calculated for the model parameters.

```
samplesSummary(samples, round = 2)

##        Mean Median St.Dev. 95%CI_low 95%CI_upp
## p      0.90   0.90    0.03      0.83      0.94
## phi[1] 0.58   0.58    0.03      0.52      0.63
## phi[2] 0.50   0.50    0.06      0.39      0.61
```

Traceplots and posterior density plots are generated using
the `samplesSummary` function from the `basicMCMCplots` package.

```
basicMCMCplots::samplesPlot(samples, legend.location = "topleft")
```

![](samplesPlot.pdf)




# Statement of need

HMC is recognized as a state-of-the-art MCMC strategy.
A testimony to this, software packages such as
`Stan` have been built exclusively around HMC sampling.
As a result, however, such software cannot operate on models
with discrete parameters where HMC cannot operate.  Models with
discrete parameters arise in a range of statistical motifs
including hidden Markov models, finite mixture models, and generally in
the presence of unobserved categorical data [@bartolucci2022discrete].
In contrast, other mainstream MCMC packages 
*`WinBUGS`, `OpenBUGS` and `jags`) can sample discrete parameters,
but provide no facilities for HMC sampling.  This leaves a gap, as there is no support
for applying HMC sampling to continuous-valued parameters of
hierarchical models which also contain discrete parameters.

`nimbleHMC` fills this gap, by providing an HMC sampler which operates
 inside `nimble`'s  MCMC engine.  `nimble` provides
a host of MCMC sampling algorithms which are suitable for either
continuous or discrete parameters, as well as the
ability to customize an MCMC algorithm by specifying sampler
assignments. `nimbleHMC` supplements the suite
of sampling algorithms provided with `nimble` with an
HMC sampler, which can be used alongside other samplers.  The example presented herein
demonstrates precisely this use case: HMC sampling operating
alongside other discrete samplers, which is not possible without the use of `nimbleHMC`.

It is an open question of what combination of
samplers will optimize MCMC efficiency.
One metric of comparison is the effective sample size of the
samples generated per unit
runtime of the algorithm.  That is, how quickly an MCMC algorithm generates
 information about the parameters.  This
metric is studied in @turek2017automated and
@ponisio2020one, but without any conclusive result.  For that reason, the ability to mix-and-match
samplers from a large pool of candidates is
important from both practical and theoretical standpoints.
Indeed, packages such as `compareMCMCs`
[@de2022comparemcmcs] are designed to compare the relative
performance of MCMC algorithms.  The addition of HMC sampling
provided by `nimbleHMC` supports new combinations of MCMC algorithms, as well as facilitates
a deeper study of practical Bayesian modelling.




# References 





