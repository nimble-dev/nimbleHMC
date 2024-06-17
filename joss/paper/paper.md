---
title: "`nimbleHMC`: An R package for Hamiltonian Monte Carlo sampling in `nimble`"
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
date: 14 January 2024
bibliography: paper.bib
output: pdf_document
---






# Summary

Markov chain Monte Carlo (MCMC) algorithms are widely used for
fitting hierarchical models to data.
MCMC is the predominant tool used in Bayesian
analyses to generate samples from the posterior
distribution of model parameters conditional on observed data.
MCMC is not a single algorithm, but rather a framework in which
various sampling methods (samplers) are assigned to operate on
subsets of unobserved parameters.
There exists a vast set of valid samplers to draw upon,
which differ in complexity, autocorrelation of samples
produced, and applicability.

Hamiltonian Monte Carlo
[HMC; @neal2011mcmc] sampling is one such technique, applicable
to continuous-valued parameters, which uses gradients to generate 
large transitions in parameter space.  The resulting samples have low
autocorrelation, and therefore have high information content,
relative for example to an equal-length sequence of highly
autocorrelated samples. The No-U-Turn (NUTS) variety of HMC sampling
[HMC-NUTS; @hoffman2014no] greatly increases the usability of HMC by
introducing a recursive tree of numerical integration steps that makes
it unnecessary to pre-specify a fixed number of steps.
@hoffman2014no also introduce a self-tuning scheme for the
step size, resulting in a fully automated HMC sampler with no need for manual tuning.

Many software packages offer implementations of MCMC, such as
`nimble` [@de2017programming], `WinBUGS` [@lunn2000winbugs], `jags` [@plummer2003jags], `pyMC`
[@fonnesbeck2015pymc], `NumPyro` [@phan2019composable], `TensorFlow
Probability` [@pang2020deep], and `Stan` [@carpenter2017stan], among others.
These packages differ, however, in their approaches to sampler
assignments.  As sampling techniques vary in
computation and quality of the samples, the
effectiveness of the MCMC algorithms will vary depending on the
software and model.

Among MCMC software packages, `nimble` is distinct in allowing easy
customization of sampler assignments from a high-level interface.
Users may assign any valid samplers to each parameter or group of parameters, selecting from samplers provided with `nimble` or samplers they have written in `nimble`'s algorithm programming system.
Samplers provided with `nimble` include random walk Metropolis-Hastings sampling [@robert1999metropolis],
slice sampling [@neal2003slice], elliptical slice sampling
[@murray2010elliptical], automated factor slice sampling [@tibbits2014automated],
conjugate sampling [@george1993conjugate], and others.

The `nimbleHMC` package provides implementations of two versions of
HMC-NUTS sampling for use within `nimble`, both written in `nimble`'s
algorithm programming system within R. Specifically, `nimbleHMC`
provides the original ("classic") HMC-NUTS algorithm as developed in
@hoffman2014no, and a modern version of
HMC-NUTS sampling matching the HMC sampler available in version 2.32.2 of `Stan`
[@stan2023stan].  The samplers provided in `nimbleHMC` can be assigned
to any continuous-valued parameters, and may be used in combination with other
samplers provided with `nimble`.





# Example

The following example demonstrates fitting a hierarchical
model to data using `nimbleHMC`.
We use the European Dipper \emph{(Cinclus cinclus)} dataset drawn from ecological
capture-recapture
[\emph{e.g.}, @lebreton1992modeling; @turek2016efficient].
Modelling includes both continuous parameters to undergo HMC
sampling and discrete parameters that cannot be sampled via HMC.
We are not aware of other software that supports this combination other than `nimbleHMC`.

Individual birds are captured, tagged, and potentially recaptured on
subsequent sighting occasions.  Data is a $255 \times 7$
binary-valued array of capture histories of 255 uniquely tagged birds
over 7 years.  Model parameters are detection
probability ($p$), and annual survival rates on non-flood years ($\phi_1$)
and flood years ($\phi_2$).  Data is provided in the R package `mra` [@mcdonald2018mra],
and individuals which are first sighted on the final (7$^{th}$) sighting occasion do not contribute 
to inference, and are removed from the sighting histories.

```
library(mra) 
data(dipper.data) 
dipper <- dipper.data[,1:7]
y <- dipper[apply(dipper, 1, which.max) < 7, ]
```

We specify the hierarchical model using uniform priors on the interval
$[0,1]$ for all parameters.
Binary-valued latent states $x_{i,t}$ represent the true alive (1) or dead (0)
state of individual $i$ on year $t$.  Doing so allows the survival
process to be modelled as $x_{i,t+1}~\sim~\text{Bernoulli}(\phi_{f_t} \cdot
x_{i,t})$ where $f_t$ indicates the flood/non-flood history of year
$t$.  The model structure conditions on the first observation of each individual, where $\text{first}_i$ is the first observation period of individual $i$, and $x_{i,\text{first}_i}$ is assigned the value one.
Observations are modelled as $y_{i,t}~\sim~\text{Bernoulli}(p \cdot x_{i,t})$.

```
library(nimbleHMC) 

code <- nimbleCode({
    phi[1] ~ dunif(0, 1)
    phi[2] ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for(i in 1:N) {
        x[i,first[i]] <- 1
        for(t in (first[i]+1):T) {
            x[i,t] ~ dbern(phi[f[t]] * x[i,t-1])
            y[i,t] ~ dbern(p * x[i,t])
        }
    }
})
```

A `nimble` model object is now built.
The argument `buildDerivs = TRUE` results in under-the-hood support for obtaining
derivatives from model calculations, as necessary for derivative-based
HMC sampling.

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
the model parameters.  `replaceSamplers` replaces the samplers operating on
$\phi_1$, $\phi_2$ and $p$ with the modern HMC-NUTS sampler (called
the `NUTS` sampler) provided in `nimbleHMC`.  The classic version of
the HMC-NUTS sampler could be assigned by specifying `type = "NUTS_classic"`.

```
conf$replaceSamplers(target = c("phi", "p"), type = "NUTS")
conf$printSamplers(byType = TRUE)

## NUTS sampler (1)
##   - phi, p 
## binary sampler (848)
##   - x[]  (848 elements)
```

Alternatively, the convenience function `configureHMC`
could be used to create an identical MCMC configuration, applying
HMC-NUTS sampling to $\phi_1$, $\phi_2$ and $p$, and default binary
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
```

Finally, posterior summary statistics are calculated for the model parameters.

```
samplesSummary(samples, round = 2)

##        Mean Median St.Dev. 95%CI_low 95%CI_upp
## p      0.90   0.90    0.03      0.84      0.95
## phi[1] 0.58   0.58    0.03      0.52      0.63
## phi[2] 0.50   0.50    0.06      0.39      0.61
```

Traceplots and posterior density plots are generated using
the `samplesPlot` function from the `basicMCMCplots` package.

```
basicMCMCplots::samplesPlot(samples, legend.location = "topleft")
```

![](samplesPlot.pdf)




# Statement of need

HMC is recognized as a state-of-the-art MCMC sampling algorithm.
As testimony to this, software packages such as
`Stan` exclusively employ HMC sampling.
Consequently, such software cannot operate on models
containing discrete parameters (upon which HMC cannot operate).  Models with
discrete parameters arise in a range of statistical motifs
including hidden Markov models, finite mixture models, and generally in
the presence of unobserved categorical data [@bartolucci2022discrete].
In contrast, other mainstream MCMC packages
(`WinBUGS`, `OpenBUGS` and `jags`) can sample discrete parameters,
but provide no facilities for HMC sampling.  This leaves a gap, as there is no support
for HMC sampling of hierarchical models that also contain discrete parameters.

`nimbleHMC` fills this gap, by providing two HMC samplers that operate
 inside `nimble`'s MCMC engine.  The base `nimble` package provides
a variety of MCMC sampling algorithms, as well as the
ability to customize MCMC sampler assignments. `nimbleHMC` augments the set
of sampling algorithms provided in `nimble` with two options for
HMC sampling, which can be used alongside any other samplers.  The example presented here
demonstrates precisely that: HMC sampling operating
alongside discrete samplers, which is not possible without the use of `nimbleHMC`.

Which combination of samplers will optimize MCMC efficiency for any particular problem is an open question.
One metric of comparison is the effective sample size of the
samples generated per unit runtime of the algorithm, which quantifies
how quickly an MCMC algorithm 
generates information about parameter posteriors.
 This metric is studied in @turek2017automated and
@ponisio2020one, with the conclusion that the best sampling strategy
is problem-specific rather than universal.  For that reason, the ability to mix-and-match
samplers from a large pool of candidates is
important from both practical and theoretical standpoints.
Indeed, packages such as `compareMCMCs`
[@de2022comparemcmcs] exist specifically to compare the relative
performance of MCMC algorithms.  The addition of HMC sampling
provided by `nimbleHMC` supports new practical combinations for applied MCMC, as well as facilitates
a deeper study of Bayesian modelling.




# References 





