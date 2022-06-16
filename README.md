# nimbleHMC

[![tests](https://github.com/nimble-dev/nimbleHMC/workflows/tests/badge.svg)](https://github.com/nimble-dev/nimbleHMC/actions)

Provides derivative-based MCMC sampling algorithms for use in conjunction with the nimble package.  These include:

- Hamiltonian Monte Carlo (HMC-NUTS) sampler
- Langevin sampler (*under development*)

See the [nimble website](https://r-nimble.org/) for more information
and examples.

<!--
The nimbleHMC package must be used with nimble version XXXX or 
higher. To check the current version number of nimble use `packageVersion("nimble")`. 
-->

The `nimbleHMC` package must be used with nimble the `"AD-rc0"` branch of
the `nimble` package.  To install this version of `nimble` and the
`nimbleHMC` package, use:

```
library(devtools)
install_github("nimble-dev/nimble", ref = "AD-rc0", subdir = "packages/nimble")
install_github("nimble-dev/nimbleHMC", subdir = "nimbleHMC")
```
