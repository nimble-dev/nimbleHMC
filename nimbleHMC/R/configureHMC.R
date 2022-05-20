#' add HMC sampler
#'
#' Add a Hamiltonian Monte Carlo sampler to an existing nimble MCMC configuration
#'
#' @param MCMCconf An existing nimble MCMC configuration, as returned by `configureMCMC`.
#' @param nodes A character vector of stochastic node names to sample by HMC.
#' If nodes is `character()`, all nodes will be sampled.
#' @param control Optional list of control parameters to be passed as the `control`
#' argument to `sampler_HMC`.  See `help(sampler_HMC)`.
#' @param replace If `TRUE`, remove any existing samplers assigned to `nodes` in `MCMCconf`
#' before adding the HMC sampler.
#'
#' @details
#'
#' This is a helper function to invoke `MCMCconf$addSampler` to add an HMC sampler to an MCMC configuration.
#'
#' Use this function if you have created an MCMC configuration and want modify it by
#' adding an MCMC sampler,
#'
#' Use `MCMCconf$addSampler` instead if you need more fine-grained control.  See `help(configureMCMC)` in nimble.
#'
#' @examples
#' nimbleOptions(enableDerivs = TRUE)
#'
#' code <- nimbleCode({
#'     b0 ~ dnorm(0, 0.001)
#'     b1 ~ dnorm(0, 0.001)
#'     sigma ~ dunif(0, 10000)
#'     for(i in 1:N) {
#'         mu[i] <- b0 + b1 * x[i]
#'         y[i] ~ dnorm(mu[i], sd = sigma)
#'     }
#' })
#'
#' N <- 10
#' constants <- list(N = N, x = 1:N)
#' data <- list(y = 1:N)
#' inits <- list(b0 = 1, b1 = 0.1, sigma = 1)
#'
#' Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
#'
#' conf <- configureMCMC(Rmodel, nodes = NULL)
#'
#' addHMC(conf) # modify conf by adding HMC sampler
#'
#' Rmcmc <- buildMCMC(conf)
#' # Cmodel <- compileNimble(Rmodel)
#' # Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
#' # samples <- runMCMC(Cmcmc)
#' @export
addHMC <- function(MCMCconf, nodes = character(), control = list(),
                   replace = FALSE) {
  model <- MCMCconf$model
  if(identical(nodes, character()))
    nodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  if(replace)
    MCMCconf$removeSamplers(nodes)
  MCMCconf$addSampler(target = nodes, type = 'HMC',  control = control)
  MCMCconf
}

#' Create HMC configuration
#'
#' Create a nimble MCMC configuration object containing only an HMC
#'
#' @param model A nimble model, such as returned by `nimbleModel`.
#' @param nodes A character vector of stochastic node names to sample by HMC.
#' If nodes is `character()`, all nodes will be sampled.
#' @param control Optional list of control parameters to be passed as the `control`.
#' argument to `sampler_HMC`.  See `help(sampler_HMC)`.
#' @param ... Other arguments that will be passed to `configureMCMC`.
#'
#' @details
#'
#' This function can be used like `configureMCMC` in nimble.  It will return
#' an MCMC configuration with only an HMC sampler.
#'
#' Use this function if you want to create an MCMC configuration and the modify it
#' before building it.  If you want to go straight to a built MCMC sampler containing
#' just an HMC sampler, use buildMCMC instead.
#'
#' @export
#'
#' @examples
#' # See example for help(addHMC)
#' # Replace the configureMCMC and addHMC lines with
#' # conf <- configureHMC(Rmodel)
configureHMC <- function(model, nodes = character(), control = list(), ...) {
  MCMCconf <- configureMCMC(model, nodes = NULL, ..., print=FALSE)
  if(is.list(control[['HMC']])) control <- control[['HMC']] # in case it is a list of lists
  addHMC(MCMCconf=MCMCconf, nodes=nodes, control=control)
  MCMCconf
}

#' Create HMC
#'
#' Create a nimble MCMC object containing only an HMC samplers.R
#' @param model A nimble model, such as returned by `nimbleModel`.
#' @param nodes A character vector of stochastic node names to sample by HMC.
#' If nodes is `character()`, all nodes will be sampled.
#' @param control Optional list of control parameters to be passed as the `control`.
#' argument to `sampler_HMC`.  See `help(sampler_HMC)`.
#' @param ... Other arguments that will be passed to `configureHMC` and from there to `configureMCMC`.
#'
#' @details
#'
#' This is the fastest way to create an HMC in nimble.  Use this if you don't want to
#' modify the MCMC configuration before or after adding an HMC sampler to it.
#'
#' @export
#' @examples
#' nimbleOptions(enableDerivs = TRUE)
#'
#' code <- nimbleCode({
#'     b0 ~ dnorm(0, 0.001)
#'     b1 ~ dnorm(0, 0.001)
#'     sigma ~ dunif(0, 10000)
#'     for(i in 1:N) {
#'         mu[i] <- b0 + b1 * x[i]
#'         y[i] ~ dnorm(mu[i], sd = sigma)
#'     }
#' })
#'
#' N <- 10
#' constants <- list(N = N, x = 1:N)
#' data <- list(y = 1:N)
#' inits <- list(b0 = 1, b1 = 0.1, sigma = 1)
#' Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
#'
#' Rmcmc <- buildHMC(Rmodel)
#'
#' # Cmodel <- compileNimble(Rmodel)
#' # Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
#' # samples <- runMCMC(Cmcmc)
buildHMC <- function(model, nodes = character(), control = list(), ...) {
  MCMCconf <- configureHMC(model=model, nodes=nodes, control = control, ...)
  buildMCMC(MCMCconf)
}
