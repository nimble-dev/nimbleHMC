

#' Add HMC sampler
#'
#' Add a No-U-Turn (NUTS) Hamiltonian Monte Carlo (HMC) sampler to an existing nimble MCMC configuration object
#'
#' @param conf A nimble MCMC configuration object, as returned by `configureMCMC`.
#' @param target A character vector of continuous-valued stochastic node names to sample.  If this argument contains any discrete-valued nodes, an error is produced and no HMC sampler is added.  If this argument is omitted, then no HMC sampler is added.
#' @param type A character string specifying the type of HMC sampler to add, either "NUTS" or "NUTS_classic".  See `help(NUTS)` or `help(NUTS_classic)` for details of each sampler.  The default sampler type is "NUTS".
#' @param control Optional named list of control parameters to be passed as the `control` argument to the HMC sampler.  See `help(NUTS)` or `help(NUTS_classic)` for details of the control list elements accepted by each sampler.
#' @param replace Logical argument.  If `TRUE`, any existing samplers operating on the specified nodes will be removed, prior to adding the HMC sampler.  Default value is `FALSE`.
#' @param print Logical argument whether to print the newly added HMC sampler.  Default value is `TRUE`.
#' @param ... Additional named arguments passed through ... will be used as additional control list elements.
#'
#' @details
#'
#' This function adds an HMC sampler to an MCMC configuration object.  Use this function if you have already created an MCMC configuration and want to add an HMC sampler.  Optionally, using `replace = TRUE`, this function will also remove any existing samplers operating on the target node(s).
#'
#' Either the `NUTS_classic` or the `NUTS` sampler can be added.  Both implement variants of No-U-Turn HMC sampling, however the `NUTS` sampler uses more modern adapatation techniques.  See `help(NUTS)` or `help(NUTS_classic)` for details.
#'
#' Use `conf$addSampler` instead if you need more fine-grained control.  See `help(configureMCMC)` in nimble.
#'
#' @author Daniel Turek
#'
#' @export
#'
#' @return Invisibly returns an object of class `MCMCconf`, but this function is primary called for its side effect.
#'
#' @seealso \code{\link{configureHMC}} \code{\link{buildHMC}} \code{\link[nimble]{configureMCMC}} \code{\link[nimble]{addSampler}} \code{\link{sampler_NUTS}} \code{\link{sampler_NUTS_classic}}
#' 
#' @examples
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
#' ## create default MCMC configuration object
#' conf <- configureMCMC(Rmodel)
#'
#' ## remove default samplers operating on b0 and b1
#' conf$removeSamplers(c("b0", "b1"))
#'
#' ## add an HMC sampler operating on b0 and b1
#' addHMC(conf, target = c("b0", "b1"))
#'
#' Rmcmc <- buildMCMC(conf)
#' 
#' # Cmodel <- compileNimble(Rmodel)
#' # Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
#' # samples <- runMCMC(Cmcmc)
addHMC <- function(conf, target = character(), type = 'NUTS', control = list(), replace = FALSE, print = TRUE, ...) {
    if(identical(target, character()))  return()
    if(!(type %in% c('NUTS', 'NUTS_classic')))   stop('type argument must be \"NUTS\" or \"NUTS_classic\".')
    targetExpanded <- conf$model$expandNodeNames(target)
    if(length(targetExpanded)) {
        discreteNodes <- targetExpanded[conf$model$isDiscrete(targetExpanded)]
        if(length(discreteNodes))   stop(paste0('HMC sampler cannot be applied to discrete nodes: ', paste0(discreteNodes, collapse = ', ')), call. = FALSE)
    }
    if(replace)   conf$removeSamplers(target)
    control <- c(control, list(...))
    conf$addSampler(target = target, type = type, control = control, print = print)
    return(invisible(conf))
}



#' Configure HMC
#'
#' Create a nimble MCMC configuration object which applies HMC sampling to continuous-valued dimensions
#'
#' @param model A nimble model, as returned by `nimbleModel`
#' @param nodes A character vector of stochastic node names to be sampled. If an empty character vector is provided (the default), then all stochastic non-data nodes will be sampled.  An HMC sampler will be applied to all continuous-valued non-data nodes, and nimble's default sampler will be assigned for all discrete-valued samplin to apply, either "NUTS" or "NUTS_classic".
#' @param type A character string specifying the type of HMC sampling to apply, either "NUTS" or "NUTS_classic".  See `help(NUTS)` or `help(NUTS_classic)` for details of each sampler.  The default sampler type is "NUTS".
#' @param control Optional named list of control parameters to be passed as the `control` argument to the HMC sampler.  See `help(NUTS)` or `help(NUTS_classic)` for details of the control list elements accepted by each sampler.
#' @param print Logical argument specifying whether to print the montiors and samplers.  Default is TRUE.
#' @param ... Other arguments that will be passed to `configureMCMC`
#'
#' @details
#'
#' This function can be used like `configureMCMC` in nimble to create an MCMC configuration object.  It will return an MCMC configuration with an HMC sampler assigned to continuous-valued model dimensions, and nimble's default sampler assigned for discrete-valued dimensions (or, only for the nodes specified in the `nodes` argument).  The resulting MCMC configuration object can be used as an argument to `buildMCMC` to generate an executable MCMC algorithm.
#'
#' Either the `NUTS_classic` or the `NUTS` sampler can be applied.  Both implement variants of No-U-Turn HMC sampling, however the `NUTS` sampler uses more modern adapatation techniques. See `help(NUTS)` or `help(NUTS_classic)` for details.
#'
#' Use this function if you want to create an MCMC configuration, and then modify it further before building the MCMC algorithm.  `buildHMC` provides a more direct route to a compilable MCMC algorithm with HMC sampling applied to all continuous-valued dimensions.
#'
#' @author Daniel Turek
#'
#' @export
#'
#' @return An object of class `MCMCconf`.
#'
#' @seealso \code{\link{addHMC}} \code{\link{buildHMC}} \code{\link[nimble]{configureMCMC}} \code{\link[nimble]{addSampler}} \code{\link{sampler_NUTS}} \code{\link{sampler_NUTS_classic}}
#' 
#' @examples
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
#' ## create MCMC configuration object with only an HMC sampler
#' conf <- configureHMC(Rmodel)
#'
#' Rmcmc <- buildMCMC(conf)
#' 
#' # Cmodel <- compileNimble(Rmodel)
#' # Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
#' # samples <- runMCMC(Cmcmc)
configureHMC <- function(model, nodes = character(), type = 'NUTS', control = list(), print = TRUE, ...) {
    if(!is.model(model)) stop('\'model\' argument must be a nimble model object',  call. = FALSE)
    if(!isTRUE(model$modelDef[['buildDerivs']])) stop('must set buildDerivs = TRUE when building model',  call. = FALSE)
    nodesProvided <- !identical(nodes, character())
    if(nodesProvided) {
        nodes <- model$expandNodeNames(nodes)
        isDiscreteBool <- model$isDiscrete(nodes)
        stochContNodes <- nodes[!isDiscreteBool]
        stochDiscNodes <- nodes[ isDiscreteBool]
    } else {
        nodeLists <- hmc_determineNodeLists(model)
        stochContNodes <- nodeLists$stochCont
        stochDiscNodes <- nodeLists$stochDisc
        postPredNodes  <- nodeLists$postPred
    }
    conf <- configureMCMC(model, nodes = NULL, print = FALSE, ...)
    addHMC(conf, stochContNodes, type, control, print = FALSE)
    conf$addSampler(target = stochDiscNodes, control = control, print = FALSE, default = TRUE)
    if(!nodesProvided)   conf$addSampler(target = postPredNodes, control = control, print = FALSE, default = TRUE)
    if(print)   conf$show()     ## conf$show(...) is needed for includeConfGetUnsampledNodes argument
    return(invisible(conf))
}


#' Build HMC
#'
#' Build an MCMC algorithm which applies HMC sampling to continuous-valued dimensions
#' 
#' @param model A nimble model, as returned by `nimbleModel`.  Alternatively, an MCMC configuration object, as returned by either `configureHMC` or `configureHMC`.  See details.
#' @param nodes A character vector of stochastic node names to be sampled. If an empty character vector is provided (the default), then all stochastic non-data nodes will be sampled.  An HMC sampler will be applied to all continuous-valued non-data nodes, and nimble's default sampler will be assigned for all discrete-valued nodes.
#' @param type A character string specifying the type of HMC sampling to apply, either "NUTS" or "NUTS_classic".  See `help(NUTS)` or `help(NUTS_classic)` for details of each sampler.  The default sampler type is "NUTS".
#' @param control Optional named list of control parameters to be passed as the `control` argument to the HMC sampler.  See `help(NUTS)` or `help(NUTS_classic)` for details of the control list elements accepted by each sampler.
#' @param print Logical argument specifying whether to print the montiors and samplers.  Default is TRUE.
#' @param ... Other arguments that will be passed to `configureHMC`.
#'
#' @details
#'
#' This is the most direct way to create an MCMC algorithm using HMC sampling in nimble.  This will create a compilable, executable MCMC algorithm, with HMC sampling assigned to all continuous-valued model dimensions, and nimble's default sampler assigned to all discrete-valued dimensions.  The `nodes` argument can be used to control which model nodes are assigned samplers.  Use this if you don't otherwise need to modify the MCMC configuration.
#'
#' Either the `NUTS_classic` or the `NUTS` sampling can be applied, which is controled by the `type` argument.  Both implement variants of No-U-Turn HMC sampling, however the `NUTS` sampler uses more modern adapatation techniques. See `help(NUTS)` or `help(NUTS_classic)` for details.
#'
#' Note that when an MCMC configuration object is provided as the `model` argument, then an executable MCMC algorithm is generated using the MCMC configuration that was provided - regardless of whether or not it specifies any HMC samplers.
#'
#' @export
#'
#' @return An object of class `MCMC`.
#'
#' @author Daniel Turek
#' 
#' @seealso \code{\link{addHMC}} \code{\link{configureHMC}} \code{\link[nimble]{configureMCMC}} \code{\link[nimble]{addSampler}} \code{\link{sampler_NUTS}} \code{\link{sampler_NUTS_classic}}
#' 
#' @examples
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
buildHMC <- function(model, nodes = character(), type = 'NUTS', control = list(), print = TRUE, ...) {
    if(is.model(model)) {
        conf <- configureHMC(model, nodes, type, control, print, includeConfGetUnsampledNodes = FALSE, ...)
    } else if(inherits(model, 'MCMCconf')) {
        conf <- model
    } else stop('\'model\' argument must be either a nimble model or an MCMC configuration object',  call. = FALSE)
    Rmcmc <- suppressMessages(buildMCMC(conf))
    return(Rmcmc)
}


#' Builds and executes NIMBLE's HMC sampler
#'
#' \code{nimbleHMC} is the most direct entry point to using NIMBLE's HMC sampler.  HMC sampling is applied to all unobserved dimensions of a hierarchical model.  Discrete-valued model dimensions cannot be sampled using HMC, and will produce an error.  See \code{help(HMC)} for details of the HMC algorithm.
#'
#' nimbleHMC provides capability for running multiple MCMC chains, specifying the number of MCMC iterations, thinning, and burn-in, and which model variables should be monitored.  It also provides options to return the posterior samples, to return summary statistics calculated from the posterior samples, and to return a WAIC value.
#'
#' The entry point for this function is providing the \code{code}, \code{constants}, \code{data} and \code{inits} arguments, to create a new NIMBLE model object, or alternatively providing an exisiting NIMBLE model object as the \code{model} argument.
#'
#' @param code The quoted code expression representing the model, such as the return value from a call to \code{nimbleCode}). Not required if \code{model} is provided.  
#' 
#' @param constants Named list of constants in the model.  Constants cannot be subsequently modified. For compatibility with JAGS and BUGS, one can include data values with constants and \code{nimbleModel} will automatically distinguish them based on what appears on the left-hand side of expressions in \code{code}.
#' 
#' @param data Named list of values for the data nodes.  Values that are NA will not be flagged as data.
#'
#' @param inits Argument to specify initial values for each MCMC chain.  See details.
#'
#' @param dimensions Named list of dimensions for variables.  Only needed for variables used with empty indices in model code that are not provided in constants or data.
#'
#' @param model A compiled or uncompiled NIMBLE model object.  When provided, this model will be used to configure the MCMC algorithm to be executed, rather than using the \code{code}, \code{constants}, \code{data} and \code{inits} arguments to create a new model object.  However, if also provided, the \code{inits} argument will still be used to initialize this model prior to running each MCMC chain.
#'
#' @param type A character string specifying the type of HMC sampling to apply, either "NUTS" or "NUTS_classic".  See `help(NUTS)` or `help(NUTS_classic)` for details of each sampler.  The default sampler type is "NUTS".
#' 
#' @param monitors A character vector giving the node names or variable names to monitor.  The samples corresponding to these nodes will returned, and/or will have summary statistics calculated. Default value is all top-level stochastic nodes of the model.
#' 
#' @param thin Thinning interval for collecting MCMC samples.  Thinning occurs after the initial nburnin samples are discarded. Default value is 1.
#' 
#' @param niter Number of MCMC iterations to run.  Default value is 10000.
#' 
#' @param nburnin Number of initial, pre-thinning, MCMC iterations to discard.  Default value is 0.
#' 
#' @param nchains Number of MCMC chains to run.  Default value is 1.
#' 
#' @param check Logical argument, specifying whether to check the model object for missing or invalid values.  Default value is \code{TRUE}.
#' 
#' @param setSeed Logical or numeric argument.  If a single numeric value is provided, R's random number seed will be set to this value at the onset of each MCMC chain.  If a numeric vector of length \code{nchains} is provided, then each element of this vector is provided as R's random number seed at the onset of the corresponding MCMC chain.  Otherwise, in the case of a logical value, if \code{TRUE}, then R's random number seed for the ith chain is set to be \code{i}, at the onset of each MCMC chain.  Note that specifying the argument \code{setSeed = 0} does not prevent setting the RNG seed, but rather sets the random number generation seed to \code{0} at the beginning of each MCMC chain.  Default value is \code{FALSE}.
#'
#' @param progressBar Logical argument.  If \code{TRUE}, an MCMC progress bar is displayed during execution of each MCMC chain.  Default value is defined by the nimble package option MCMCprogressBar.
#'
#' @param samples Logical argument.  If \code{TRUE}, then posterior samples are returned from each MCMC chain.  These samples are optionally returned as \code{coda} \code{mcmc} objects, depending on the \code{samplesAsCodaMCMC} argument.  Default value is \code{TRUE}.  See details.
#'
#' @param samplesAsCodaMCMC Logical argument.  If \code{TRUE}, then a \code{coda} \code{mcmc} object is returned instead of an R matrix of samples, or when \code{nchains > 1} a \code{coda} \code{mcmc.list} object is returned containing \code{nchains} \code{mcmc} objects.  This argument is only used when \code{samples} is \code{TRUE}.  Default value is \code{FALSE}.  See details.
#' 
#' @param summary Logical argument.  When \code{TRUE}, summary statistics for the posterior samples of each parameter are also returned, for each MCMC chain.  This may be returned in addition to the posterior samples themselves.  Default value is \code{FALSE}.  See details.
#' 
#' @param WAIC Logical argument.  When \code{TRUE}, the WAIC (Watanabe, 2010) of the model is calculated and returned.  If multiple chains are run, then a single WAIC value is calculated using the posterior samples from all chains.  Default value is \code{FALSE}. Note that the version of WAIC used is the default WAIC conditional on random effects/latent states and without any grouping of data nodes. See \code{help(waic)} for more details.
#' 
#' @param userEnv Environment in which if-then-else statements in model code will be evaluated if needed information not found in \code{constants}; intended primarily for internal use only.
#' 
#' @return A list is returned with named elements depending on the arguments, unless only one among samples, summary, and WAIC are requested, in which case only that element is returned.  These elements may include \code{samples}, \code{summary}, and \code{WAIC}.  When \code{nchains = 1}, posterior samples are returned as a single matrix, and summary statistics as a single matrix.  When \code{nchains > 1}, posterior samples are returned as a list of matrices, one matrix for each chain, and summary statistics are returned as a list containing \code{nchains+1} matrices: one matrix corresponding to each chain, and the final element providing a summary of all chains, combined.  If \code{samplesAsCodaMCMC} is \code{TRUE}, then posterior samples are provided as \code{coda} \code{mcmc} and \code{mcmc.list} objects.  When \code{WAIC} is \code{TRUE}, a WAIC summary object is returned.
#'
#' @details
#'
#' At least one of \code{samples}, \code{summary} or \code{WAIC} must be \code{TRUE}, since otherwise, nothing will be returned.  Any combination of these may be \code{TRUE}, including possibly all three, in which case posterior samples, summary statistics, and WAIC values are returned for each MCMC chain.
#'
#' When \code{samples = TRUE}, the form of the posterior samples is determined by the \code{samplesAsCodaMCMC} argument, as either matrices of posterior samples, or \code{coda} \code{mcmc} and \code{mcmc.list} objects.
#'
#' Posterior summary statistics are returned individually for each chain, and also as calculated from all chains combined (when \code{nchains > 1}).
#'
#' The \code{inits} argument can be one of three things:
#' 
#' (1) a function to generate initial values, which will be executed once to initialize the model object, and once to generate initial values at the beginning of each MCMC chain, or
#' (2) a single named list of initial values which, will be used to initialize the model object and for each MCMC chain, or
#' (3) a list of length \code{nchains}, each element being a named list of initial values.  The first element will be used to initialize the model object, and once element of the list will be used for each MCMC chain.
#' 
#' The \code{inits} argument may also be omitted, in which case the model will not be provided with initial values.  This is not recommended.
#'
#' The \code{niter} argument specifies the number of pre-thinning MCMC iterations, and the \code{nburnin} argument specifies the number of pre-thinning MCMC samples to discard.  After discarding these burn-in samples, thinning of the remaining samples will take place.  The total number of posterior samples returned will be floor((niter-nburnin)/thin).
#' 
#' @examples
#' 
#' \donttest{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, sd = 1000)
#'     sigma ~ dunif(0, 1000)
#'     for(i in 1:10) {
#'         x[i] ~ dnorm(mu, sd = sigma)
#'     }
#' })
#' data <- list(x = c(2, 5, 3, 4, 1, 0, 1, 3, 5, 3))
#' inits <- function() list(mu = rnorm(1,0,1), sigma = runif(1,0,10))
#' mcmc.output <- nimbleHMC(code, data = data, inits = inits,
#'                          monitors = c("mu", "sigma"), thin = 10,
#'                          niter = 20000, nburnin = 1000, nchains = 3,
#'                          summary = TRUE, WAIC = TRUE)
#' }
#'
#' @seealso \code{\link{configureHMC}} \code{\link{buildHMC}} \code{\link[nimble]{configureMCMC}} \code{\link[nimble]{buildMCMC}} \code{\link[nimble]{runMCMC}}
#' 
#' @author Daniel Turek
#'
#' @export
nimbleHMC <- function(code,
                      constants = list(),
                      data = list(),
                      inits,
                      dimensions = list(),
                      model,
                      type = 'NUTS',
                      monitors,
                      thin = 1,
                      niter = 10000,
                      nburnin = 0,
                      nchains = 1,
                      check = TRUE,
                      setSeed = FALSE,
                      progressBar = getNimbleOption('MCMCprogressBar'),
                      samples = TRUE,
                      samplesAsCodaMCMC = FALSE,
                      summary = FALSE,
                      WAIC = FALSE,
                      userEnv = parent.frame()) {
    if(missing(code) && missing(model)) stop('must provide either code or model argument')
    if(!samples && !summary && !WAIC) stop('no output specified, use samples = TRUE, summary = TRUE, or WAIC = TRUE')
    if(!missing(code) && inherits(code, 'modelBaseClass')) model <- code   ## let's handle it, if model object is provided as un-named first argument
    Rmodel <- mcmc_createModelObject(model, inits, nchains, setSeed, code, constants, data, dimensions, check, buildDerivs = TRUE, userEnv = userEnv)
    Rmcmc <- buildHMC(Rmodel, type = type, monitors = monitors, thin = thin, enableWAIC = WAIC, print = FALSE)
    compiledList <- compileNimble(Rmodel, Rmcmc)    ## only one compileNimble() call
    Cmcmc <- compiledList$Rmcmc
    runMCMC(Cmcmc, niter = niter, nburnin = nburnin, nchains = nchains, inits = inits,
            setSeed = setSeed, progressBar = progressBar, samples = samples,
            samplesAsCodaMCMC = samplesAsCodaMCMC, summary = summary, WAIC = WAIC)
}


## create the lists of model nodes for use in HMC configuration functions
hmc_determineNodeLists <- function(model) {
    ppNodes <- model$getNodeNames(predictiveOnly = TRUE)   ## stochastic only, already
    stochNodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE, includePredictive = FALSE)
    if(length(stochNodes) > 0) {
        isDiscreteBool <- model$isDiscrete(stochNodes)
        stochContNodes <- stochNodes[!isDiscreteBool]
        stochDiscNodes <- stochNodes[ isDiscreteBool]
    } else {
        stochContNodes <- character()
        stochDiscNodes <- character()
    }
    return(list(postPred  = ppNodes,
                stochCont = stochContNodes,
                stochDisc = stochDiscNodes))
}

