



#'
#' myTempSampler name
#'
#' myTempSampler description of what it does
#'
#' @import nimble
#' 
#' @param model A nimble model object.  See details.
#' @param mvSaved A nimble modelValues object.  See details.
#' @param target A character vector of node names.  See details.
#' @param control A named list of sampler control elements.  See details.
#'
#' @details XXXX DESCRIBE ARGUMENTS AND SAMPLER USE
#'
#' @author Daniel Turek
#'
#' @examples
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
#' constants <- list(N = 10)
#' Rmodel <- nimbleModel(code, constants)
#' conf <- configureMCMC(Rmodel)
#' conf$setSamplers()
#' conf$addSampler(target = c('b0', 'b1', 'sigma'), type = 'RW_block')  ## XXXXXX MAKE 'HMC'
#' Rmcmc <- buildMCMC(conf)
#' 
#' @export
myTempSampler <- nimbleFunction(
    name = 'myTempSampler',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) { },
    run = function() {
        model[[target]] <<- model[[target]] + 1
        calculate(model, target)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
    },
    methods = list(
        reset = function() { }
    )
)


