# nimbleHMC

[![tests](https://github.com/nimble-dev/nimbleHMC/workflows/tests/badge.svg)](https://github.com/nimble-dev/nimbleHMC/actions)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.06745/status.svg)](https://doi.org/10.21105/joss.06745)

Provides derivative-based MCMC sampling algorithms and convenence functions, for use in conjunction with the MCMC engine avaialble in the `nimble` package.  Sampling algorithms include:

- No-U-Turn Hamiltonian Monte Carlo (`NUTS`) sampler
- Historical implementation of the original No-U-Turn HMC (`NUTS_classic`) sampler
- Langevin sampler (*under development*)

See the HMC section of the [Nimble User Manual](https://r-nimble.org/html_manual/cha-mcmc.html#subsec:HMC) for more information and examples.

[General package information](https://cran.r-project.org/web/packages/nimbleHMC/) about `nimbleHMC`, and the [complete API for package functions](https://cran.r-project.org/web/packages/nimbleHMC/nimbleHMC.pdf) are available on CRAN.

Additional information about the `nimble`package itself is available at the [nimble website](https://r-nimble.org/).



### Installation and Package Requirements

Use of `nimbleHMC` requires installation of the core `nimble` package.  Detailed instructions for installing `nimble` are available in the [`nimble` package README](https://github.com/nimble-dev/nimble/blob/devel/README.md).

`nimbleHMC` must be used with version `1.0.0` or higher of `nimble`, or the latest version available on CRAN.  To check the version number of the currently installed version of `nimble`, use:

```r
packageVersion("nimble")
```

The `nimbleHMC` package itself can be installed directly from CRAN, using:

```r
install.packages("nimbleHMC")
```


<!--
library(remotes)
remotes::install_github("nimble-dev/nimble", ref = "devel", subdir = "packages/nimble")
remotes::install_github("nimble-dev/nimbleHMC", ref="master", subdir = "nimbleHMC")

For errors during installation of `nimbleHMC` occuring on Windows machines, relating to either of the following error messages:

Error: package 'nimble' is not installed for 'arch = i386'
Error: loading failed for 'i386'

try installing the `nimbleHMC` package using:

remotes::install_github("nimble-dev/nimbleHMC", ref="master", subdir = "nimbleHMC", INSTALL_opts=c("--no-multiarch"))
-->



### Automatic Differentiation

`nimbleHMC` makes use of the automatic differentiation (AD) feature of `nimble`, which was released in `nimble` version 1.0.0.  See [Chapter 16: Automatic Derivatives](https://r-nimble.org/html_manual/cha-AD.html) of the Nimble User Manual for more information about the capabilities of the AD system, and how to use the AD system to calculate derivatives of functions written as `nimble` algorithms, as well as derivatives of model calculations.

<!--
In order to use HMC sampling (and other derivative-based algorithms), derivatives need to be enabled for `nimble` using the setting:
nimbleOptions(enableDerivs = TRUE)
-->

For using HMC sampling on a model, derivative calculations need to be built into for the model object.  This is accomplished using the `buildDerivs = TRUE` argument in the call to `nimbleModel` as:
```
nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
```

For using HMC sampling on models that include user-defined distributions, you will need to include the argument `buildDerivs = TRUE` in the definition of your distribution, as:

```
dmy_distribution <- nimbleFunction(
    run = function(x = double(), param = double(), log = integer(0, default = 0)) {
        ...
    },
    buildDerivs = TRUE
)
```


### Contributing and Support

Contributions to the `nimbleHMC` package should be submitted via pull request on GitHub.  For additional guideliens on making contributions, please see the [contributing guidelines for the `nimble` package](https://github.com/nimble-dev/nimble/blob/devel/CONTRIBUTING.md).

Issues, feature requests, or bugs should be reported using GitHub issues, submitted to the `nimbleHMC` repository.

For additional support using `nimble` or `nimbleHMC`, please see the [Nimble User Manual](https://r-nimble.org/html_manual/cha-welcome-nimble.html).  Any additional questions can be submitted to the [nimble-users Google group](https://groups.google.com/g/nimble-users).

