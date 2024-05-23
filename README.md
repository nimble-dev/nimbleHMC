# nimbleHMC

[![tests](https://github.com/nimble-dev/nimbleHMC/workflows/tests/badge.svg)](https://github.com/nimble-dev/nimbleHMC/actions)

Provides derivative-based MCMC sampling algorithms and convenence functions, for use in conjunction with the MCMC engine avaialble in the `nimble` package.  Sampling algorithms include:

- No-U-Turn Hamiltonian Monte Carlo (`NUTS`) sampler
- Historical implementation of the original No-U-Turn HMC (`NUTS_classic`) sampler
- Langevin sampler (*under development*)

See the HMC section of the [nimble user manual](https://r-nimble.org/html_manual/cha-mcmc.html#subsec:HMC) for more information and examples.

[General package information](https://cran.r-project.org/web/packages/nimbleHMC/) about `nimbleHMC`, and the [complete API for package functions](https://cran.r-project.org/web/packages/nimbleHMC/nimbleHMC.pdf) are available on CRAN.

Additional information about the `nimble`package itself is available at the [nimble website](https://r-nimble.org/).

<!--
The nimbleHMC package must be used with nimble version XXXX or 
higher. To check the current version number of nimble use `packageVersion("nimble")`. 
-->

### Package Requirements

`nimbleHMC` must be used with version `1.0.0` or higher of the `nimble` package.

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

`nimbleHMC` makes use of the automatic differentiation (AD) feature of `nimble`, which is currently available as a beta release.  See [nimble AD beta release](https://r-nimble.org/ad-beta) for more information about models and algorithms that make use of the AD features of `nimble`.

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

For additional support using `nimble` or `nimbleHMC`, please see the [nimble user manual](https://r-nimble.org/html_manual/cha-welcome-nimble.html).  Any additional questions can be submitted to the [nimble-users Google group](https://groups.google.com/g/nimble-users).

