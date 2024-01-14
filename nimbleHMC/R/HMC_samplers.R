


## #' Langevin Sampler
## #'
## #' The langevin sampler implements a special case of Hamiltonian Monte Carlo (HMC) sampling where only a single leapfrog step is taken on each sampling iteration, and the leapfrog step-size is adapted to match the scale of the posterior distribution (independently for each dimension being sampled). The single leapfrog step is done by introducing auxiliary momentum variables, and using first-order derivatives to simulate Hamiltonian dynamics on this augmented paramter space (Neal, 2011). Langevin sampling can operate on one or more continuous-valued posterior dimensions. This sampling technique is also known as Langevin Monte Carlo (LMC), and the Metropolis-Adjusted Langevin Algorithm (MALA).
## #'
## #' @param model An uncompiled nimble model object on which the MCMC will operate.
## #' @param mvSaved A nimble \code{modelValues} object to be used to store MCMC samples.
## #' @param target A character vector of node names on which the sampler will operate.
## #' @param control A named list that controls the precise behavior of the sampler. The default values for control list elements are specified in the setup code of the sampler. A description of the possible control list elements appear in the details section.
## #'
## #' @details
## #'
## #' The Langevin sampler accepts the following control list elements:
## #'
## #' \itemize{
## #' \item scale. An optional multiplier, to scale the step-size of the leapfrog steps. If adaptation is turned off, this uniquely determines the leapfrog step-size (default = 1)
## #' \item adaptive. A logical argument, specifying whether the sampler will adapt the leapfrog step-size (scale) throughout the course of MCMC execution. The scale is adapted independently for each dimension being sampled. (default = TRUE)
## #' \item adaptInterval. The interval on which to perform adaptation. (default = 200)
## #' }
## #' 
## #' @import nimble
## #' @import methods
## #' 
## #' @aliases Langevin langevin
## #' 
## #' @author Daniel Turek
## #' 
## #' @examples
## #' code <- nimbleCode({
## #'     b0 ~ dnorm(0, 0.001)
## #'     b1 ~ dnorm(0, 0.001)
## #'     sigma ~ dunif(0, 10000)
## #'     for(i in 1:N) {
## #'         mu[i] <- b0 + b1 * x[i]
## #'         y[i] ~ dnorm(mu[i], sd = sigma)
## #'     }
## #' })
## #' 
## #' N <- 10
## #' constants <- list(N = N, x = 1:N)
## #' data <- list(y = 1:N)
## #' inits <- list(b0 = 1, b1 = 0.1, sigma = 1)
## #' 
## #' Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
## #' 
## #' conf <- configureMCMC(Rmodel, nodes = NULL)
## #' 
## #' conf$addSampler(target = c('b0', 'b1', 'sigma'), type = 'HMC')
## #' 
## #' Rmcmc <- buildMCMC(conf)

sampler_langevin <- nimbleFunction(
    name = 'sampler_langevin',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        stop('The Langevin sampler is not yet fully implemented.')  ## XXXXXXXXXXXXXX
        ## control list extraction
        scale         <- if(!is.null(control$scale))         control$scale         else 1      ## step-size multiplier
        adaptive      <- if(!is.null(control$adaptive))      control$adaptive      else TRUE
        adaptInterval <- if(!is.null(control$adaptInterval)) control$adaptInterval else 200
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes <- model$getDependencies(target)
        ## numeric value generation
        d <- length(targetAsScalar)
        scaleVec <- matrix(1, nrow = d, ncol = 1)
        epsilonVec <- scale * scaleVec
        q <- matrix(0, nrow = d, ncol = 1)
        p <- matrix(0, nrow = d, ncol = 1)
        grad <- matrix(0, nrow = d, ncol = 1)
        empirSamp <- matrix(0, nrow = adaptInterval, ncol = d)
        timesRan <- 0
        timesAdapted <- 0
        ## checks
        if(!isTRUE(nimbleOptions('enableDerivs')))   stop('must enable NIMBLE derivatives, set nimbleOptions(enableDerivs = TRUE)', call. = FALSE)
        if(!isTRUE(model$modelDef[['buildDerivs']])) stop('must set buildDerivs = TRUE when building model',  call. = FALSE)
        if(any(model$isDiscrete(targetAsScalar)))    stop(paste0('langevin sampler can only operate on continuous-valued nodes:', paste0(targetAsScalar[model$isDiscrete(targetAsScalar)], collapse=', ')), call. = FALSE)
    },
    run = function() {
        q[1:d, 1] <<- values(model, target)         ## current position variables
        for(i in 1:d)   p[i, 1] <<- rnorm(1, 0, 1)  ## randomly draw momentum variables
        currentH <- model$getLogProb(calcNodes) - sum(p^2)/2
        p <<- p + epsilonVec * jacobian(q) / 2      ## initial half step for p
        q <<- q + epsilonVec * p                    ## full step for q
        p <<- p + epsilonVec * jacobian(q) / 2      ## final half step for p
        propH <- model$calculate(calcNodes) - sum(p^2)/2
        jump <- decide(propH - currentH)
        if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        if(adaptive)     adaptiveProcedure()
    },
    methods = list(
        jacobian = function(q = double(2)) {
            values(model, target) <<- q[1:d, 1]
            derivsOutput <- nimDerivs(model$calculate(calcNodes), order = 1, wrt = target)
            grad[1:d, 1] <<- derivsOutput$jacobian[1, 1:d]
            returnType(double(2))
            return(grad)
        },
        adaptiveProcedure = function() {
            ## adapts leapfrog step-size for each dimension (epsilonVec),
            ## as described in Neal, Handbook of MCMC (2011), Chapter 5 (section 5.4.2.4)
            timesRan <<- timesRan + 1
            empirSamp[timesRan, 1:d] <<- values(model, target)
            if(timesRan %% adaptInterval == 0) {
                timesRan <<- 0
                timesAdapted <<- timesAdapted + 1
                gamma1 <- 1/((timesAdapted + 3)^0.8)
                for(i in 1:d)     scaleVec[i, 1] <<- gamma1 * sd(empirSamp[, i]) + (1-gamma1) * scaleVec[i, 1]
                epsilonVec <<- scale * scaleVec
            }
        },
        reset = function() {
            timesRan     <<- 0
            timesAdapted <<- 0
            scaleVec     <<- matrix(1, nrow = d, ncol = 1)
            epsilonVec   <<- scale * scaleVec
        }
    )
)


hmc_checkWarmup <- function(warmupMode, warmup, samplerName) {
    if(!(warmupMode %in% c('default', 'burnin', 'fraction', 'iterations')))   stop('`warmupMode` control argument of ', samplerName, ' sampler must have value "default", "burnin", "fraction", or "iterations".  The value provided was: ', warmupMode, '.', call. = FALSE)
    if(warmupMode == 'fraction')
        if(!is.numeric(warmup) | warmup < 0 | warmup > 1)   stop('When the `warmupMode` control argument of ', samplerName, ' sampler is "fraction", the `warmup` control argument must be a number between 0 and 1, which will specify the fraction of the total MCMC iterations to use as warmup.  The value provided for the `warmup` control argument was: ', warmup, '.', call. = FALSE)
    if(warmupMode == 'iterations')
        if(!is.numeric(warmup) | warmup < 0 | floor(warmup) != warmup)   stop('When the `warmupMode` control argument of ', samplerName, ' sampler is "iterations", the `warmup` control argument must be a non-negative integer, which will specify the number MCMC iterations to use as warmup.  The value provided for the `warmup` control argument was: ', warmup, '.', call. = FALSE)
}



hmc_setWarmup <- nimbleFunction(
    setup = function(warmupMode, warmup, messages, samplerName, targetNodesToPrint) {},
    run = function(MCMCniter = double(), MCMCnburnin = double(), adaptive = logical()) {
        ##
        ## set nwarmup
        if(warmupMode == 'default') {
            if(MCMCnburnin > 0)   nwarmup <- MCMCnburnin
            else                  nwarmup <- floor(MCMCniter/2)
        }
        if(warmupMode == 'burnin')       nwarmup <- MCMCnburnin
        if(warmupMode == 'fraction')     nwarmup <- floor(warmup*MCMCniter)
        if(warmupMode == 'iterations')   nwarmup <- warmup
        ##
        ## informative message
        if(messages) {
            if(!adaptive) {   ## adaptive = FALSE
                print('  [Note] ', samplerName, ' sampler (nodes: ', targetNodesToPrint, ') has adaptation turned off,\n         so no warmup period will be used.')
            } else {          ## adaptive = TRUE
                if(warmupMode == 'default') {
                    if(MCMCnburnin > 0)          print("  [Note] ", samplerName, " sampler (nodes: ", targetNodesToPrint, ") is using ", nwarmup, " warmup iterations.\n         Since `warmupMode` is 'default' and `nburnin` > 0,\n         the number of warmup iterations is equal to `nburnin`.\n         The burnin samples will be discarded, and all samples returned will be post-warmup.")
                    else                         print("  [Note] ", samplerName, " sampler (nodes: ", targetNodesToPrint, ") is using ", nwarmup, " warmup iterations.\n         Since `warmupMode` is 'default' and `nburnin` = 0,\n         the number of warmup iterations is equal to `niter/2`.\n         No samples will be discarded, so the first half of the samples returned\n         are from the warmup period, and the second half of the samples are post-warmup.")
                }
                if(warmupMode == 'burnin')
                    if(MCMCnburnin > 0)           print("  [Note] ", samplerName, " sampler (nodes: ", targetNodesToPrint, ") is using ", nwarmup, " warmup iterations.\n         Since `warmupMode` is 'burnin', the number of warmup iterations is equal to `nburnin`.\n         The burnin samples will be discarded, and all samples returned will be post-warmup.")
                    else
                        print("  [Note] ", samplerName, " sampler (nodes: ", targetNodesToPrint, ") is using 0 warmup iterations.\n         No adaptation is being done, apart from initialization of epsilon\n         (if `initializeEpsilon` is TRUE).")
                if(warmupMode == 'fraction') {
                    if(MCMCnburnin < nwarmup)    print("  [Note] ", samplerName, " sampler (nodes: ", targetNodesToPrint, ") is using ", nwarmup, " warmup iterations.\n         Since `warmupMode` is 'fraction', the number of warmup iterations is equal to\n         `niter*fraction`, where `fraction` is the value of the `warmup` control argument.\n         Because `nburnin` is less than the number of warmup iterations,\n         some of the samples returned will be collected during the warmup period,\n         and the remainder of the samples returned will be post-warmup.")
                    else                         print("  [Note] ", samplerName, " sampler (nodes: ", targetNodesToPrint, ") is using ", nwarmup, " warmup iterations.\n         Since `warmupMode` is 'fraction', the number of warmup iterations is equal to\n          `niter*fraction`, where `fraction` is the value of the warmup `control` argument.\n         Because `nburnin` exceeds the number of warmup iterations,\n         all samples returned will be post-warmup.")
                }
                if(warmupMode == 'iterations')
                    if(MCMCnburnin < nwarmup)    print("  [Note] ", samplerName, " sampler (nodes: ", targetNodesToPrint, ") is using ", nwarmup, " warmup iterations.\n         Since `warmupMode` is 'iterations', the number of warmup iterations\n         is the value of the `warmup` control argument.\n         Because `nburnin` is less than the number of warmup iterations,\n         some of the samples returned will be collected during the warmup period,\n         and the remainder of the samples returned will be post-warmup.")
                    else                         print("  [Note] ", samplerName, " sampler (nodes: ", targetNodesToPrint, ") is using ", nwarmup, " warmup iterations.\n         Since `warmupMode` is 'iterations', the number of warmup iterations\n         is the value of the `warmup` control argument.\n         Because `nburnin` exceeds the number of warmup iterations,\n         all samples returned will be post-warmup.")
            }
        }
        ##
        ## hard check that nwarmup >= 20
        if(adaptive & nwarmup > 0 & nwarmup < 20) {
            print('  [Error] ', samplerName, ' sampler (nodes: ', targetNodesToPrint, ') requires a minimum of 20 warmup iterations.')
            stop()
        }
        ##
        returnType(double())
        return(nwarmup)
    }
)



#' Classic No-U-Turn (NUTS_classic) Hamiltonian Monte Carlo (HMC) Sampler
#'
#' The NUTS_classic sampler implements the original No-U-Turn (NUTS classic) sampler as put forth in Hoffman and Gelman (2014) for performing joint updates of multiple continuous-valued posterior dimensions. This is done by introducing auxiliary momentum variables and using first-order derivatives to simulate Hamiltonian dynamics on this augmented paramter space. Internally, any posterior dimensions with bounded support are transformed, so sampling takes place on an unconstrained space. In contrast to standard HMC (Neal, 2011), the NUTS_classic algorithm removes the tuning parameters of the leapfrog step size and the number of leapfrog steps, thus providing a sampling algorithm that can be used without hand tuning or trial runs.
#'
#' @param model An uncompiled nimble model object on which the MCMC will operate.
#' @param mvSaved A nimble \code{modelValues} object to be used to store MCMC samples.
#' @param target A character vector of node names on which the sampler will operate.
#' @param control A named list that controls the precise behavior of the sampler. The default values for control list elements are specified in the setup code of the sampler. A description of the possible control list elements appear in the details section.
#'
#' @details
#'
#' The NUTS_classic sampler accepts the following control list elements:
#' 
#' \itemize{
#' \item messages. A logical argument, specifying whether to print informative messages. (default = TRUE)
#' \item numWarnings. A numeric argument, specifying how many warnings messages to emit (for example, when \code{NaN} values are encountered). See additional details below. (default = 0)
#' \item epsilon. A positive numeric argument, specifying the initial step-size value. If not provided, an appropriate initial value is selected.
#' \item gamma. A positive numeric argument, specifying the degree of shrinkage used during the initial period of step-size adaptation. (default = 0.05)
#' \item t0. A non-negative numeric argument, where larger values stabilize (attenuate) the initial period of step-size adaptation. (default = 10)
#' \item kappa. A numeric argument between zero and one, where smaller values give a higher weighting to more recent iterations during the initial period of step-size adaptation. (default = 0.75)
#' \item delta. A numeric argument, specifying the target acceptance probability used during the initial period of step-size adaptation. (default = 0.65)
#' \item deltaMax. A positive numeric argument, specifying the maximum allowable divergence from the Hamiltonian value. Paths which exceed this value are considered divergent and will not proceed further. (default = 1000)
#' \item M. A vector of positive real numbers, with length equal to the number of dimensions being sampled. Elements of \code{M} specify the diagonal elements of the diagonal mass matrix (or the metric) used for the auxiliary momentum variables in sampling. Sampling may be improved if the elements of \code{M} approximate the marginal inverse variance (precision) of the (potentially transformed) parameters. (default: a vector of ones).
#' \item warmupMode. A character string, specifying the behavior for choosing the number of warmup iterations. Four values are possible. The value 'default' (the default) sets the number of warmup iterations as the number of burnin iterations (if a positive value for \code{nburnin} is used) or half the number of MCMC iterations in each chain (if \code{nburnin = 0}). The value 'burnin' sets the number of warmup iterations as the number of burnin iterations regardless of the length of the burnin period. The value 'fraction' sets the number of warmup iterations as \code{fraction*niter}, where \code{fraction} is the value of the \code{warmup} control argument, and \code{niter} is the number of MCMC iterations in each chain; in this case, the value of the \code{warmup} control argument must be between 0 and 1. The value 'iterations' sets the number of warmup iterations as the value of the \code{warmup} control argumnet, regardless of the length of the burnin period or the number of MCMC iterations; in this case the value of \code{warmup} must be a non-negative integer. In all cases, the number of (pre-thinning) samples discarded equals \code{nburnin}, as is always the case for MCMC in NIMBLE.
#' \item warmup. Numeric value used in determining the number of warmup iterations. This control argument is only used when \code{warmupMode} is 'fraction' or 'iterations'. 
#' \item maxTreeDepth. The maximum allowable depth of the binary leapfrog search tree for generating candidate transitions. (default = 10)
#' \item adaptWindow. Number of iterations in the first adaptation window used for adapting the mass matrix (M). Subsequent adaptation windows double in length, so long as enough warmup iterations are available. (default = 25)
#' \item initBuffer. Number of iterations in the initial warmup window, which occurs prior to the first adaptation of the metric M. (default = 75)
#' \item termBuffer. Number of iterations in the final (terminal) warmup window, before which the metric M is not adjusted. (default = 50)
#' \item adaptive. A logical argument, specifying whether to do any adaptation whatsoever. When \code{TRUE}, specific adaptation routines are controlled by the \code{adaptEpsilon} and \code{adaptM} control list elements. (default = TRUE)
#' \item adaptEpsilon. A logical argument, specifying whether to perform stepsize adaptation. Only used when \code{adaptive = TRUE}. (default = TRUE)
#' \item adaptM. A logical argument, specifying whether to perform adaptation of the mass matrix (metric) M. Only used when \code{adaptive = TRUE}. (default = TRUE)
#' \item initializeEpsilon. A logical argument, specifying whether to perform the epsilon (stepsize) initialization routine at the onset of each adaptation window. (default = TRUE)
#' }
#'
#' \code{NaN} values may be encountered in the course of the leapfrog procedure. In particular, when the stepsize (\code{epsilon}) is too large, the leapfrog procedure can step too far and arrive at an invalid region of parameter space, thus generating a \code{NaN} value in the likelihood evaluation or in the gradient calculation. These situation are handled by the sampler by rejecting the \code{NaN} value, and reducing the stepsize.
#' 
#' @import nimble
#' 
#' @export
#'
#' @return A object of class `sampler_NUTS_classic`.
#' 
#' @aliases NUTS-classic NUTS_classic nuts-classic nuts_classic sampler_NUTS_classic
#' 
#' @author Daniel Turek
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
#' conf <- configureMCMC(Rmodel, nodes = NULL)
#' 
#' conf$addSampler(target = c('b0', 'b1', 'sigma'), type = 'NUTS_classic')
#' 
#' Rmcmc <- buildMCMC(conf)
#'
#' @references
#'
#' Hoffman, Matthew D., and Gelman, Andrew (2014). The No-U-Turn Sampler: Adaptively setting path lengths in Hamiltonian Monte Carlo. \emph{Journal of Machine Learning Research}, 15(1): 1593-1623.
sampler_NUTS_classic <- nimbleFunction(
    name = 'sampler_NUTS_classic',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        printTimesRan  <- extractControlElement(control, 'printTimesRan',  FALSE)
        printEpsilon   <- extractControlElement(control, 'printEpsilon',   FALSE)
        printJ         <- extractControlElement(control, 'printJ',         FALSE)
        printM         <- extractControlElement(control, 'printM',         FALSE)
        messages       <- extractControlElement(control, 'messages',       getNimbleOption('verbose'))
        numWarnings    <- extractControlElement(control, 'numWarnings',    0)
        epsilon        <- extractControlElement(control, 'epsilon',        0)          # initial epsilon, if 0 then use 1
        initializeEpsilon <- extractControlElement(control, 'initializeEpsilon', TRUE) # use initializeEpsilon step?
        gamma          <- extractControlElement(control, 'gamma',          0.05)
        t0             <- extractControlElement(control, 't0',             10)
        kappa          <- extractControlElement(control, 'kappa',          0.75)
        delta          <- extractControlElement(control, 'delta',          0.65)
        deltaMax       <- extractControlElement(control, 'deltaMax',       1000)
        M              <- extractControlElement(control, 'M',              -1)
        warmupMode     <- extractControlElement(control, 'warmupMode',     'default')   ## 'default', 'burnin', 'fraction', or 'iterations'
        warmup         <- extractControlElement(control, 'warmup',         -1)          ## used if warmupMode is 'fraction' or 'iterations'
        maxTreeDepth   <- extractControlElement(control, 'maxTreeDepth',   10)
        adaptWindow    <- extractControlElement(control, 'adaptWindow',    25)
        initBuffer     <- extractControlElement(control, 'initBuffer',     75)
        termBuffer     <- extractControlElement(control, 'termBuffer',     50)
        adaptive       <- extractControlElement(control, 'adaptive',       TRUE) # any adaptation? (if FALSE, next two flags are ignored)
        adaptEpsilon   <- extractControlElement(control, 'adaptEpsilon',   TRUE) # stepsize adaptation?
        adaptM         <- extractControlElement(control, 'adaptM',         TRUE) # mass matrix adaptation?
        ## node list generation
        targetNodes <- model$expandNodeNames(target)
        if(length(targetNodes) <= 0) stop('NUTS_classic sampler must operate on at least one node', call. = FALSE)
        targetNodesAsScalars <- model$expandNodeNames(targetNodes, returnScalarComponents = TRUE)
        targetNodesToPrint <- paste(targetNodes, collapse = ', ')
        if(nchar(targetNodesToPrint) > 100)   targetNodesToPrint <- paste0(substr(targetNodesToPrint, 1, 97), '...')
        calcNodes <- model$getDependencies(targetNodes)
        ## check for discrete nodes (early, before parameterTransform is specialized)
        if(any(model$isDiscrete(targetNodesAsScalars))) stop(paste0('NUTS_classic sampler cannot operate on discrete-valued nodes: ', paste0(targetNodesAsScalars[model$isDiscrete(targetNodesAsScalars)], collapse = ', ')))
        ## processing of bounds and transformations
        my_parameterTransform <- parameterTransform(model, targetNodesAsScalars)
        d <- my_parameterTransform$getTransformedLength()
        d2 <- max(d, 2) ## for pre-allocating vectors
        nimDerivs_wrt <- 1:d
        derivsInfo_return <- makeModelDerivsInfo(model, targetNodes, calcNodes)
        nimDerivs_updateNodes   <- derivsInfo_return$updateNodes
        nimDerivs_constantNodes <- derivsInfo_return$constantNodes
        ## numeric value generation
        timesRan <- 0;   nwarmup <- 0;   mu <- 0;   logEpsilonBar <- 0;   Hbar <- 0
        q <- numeric(d2);   qL <- numeric(d2);   qR <- numeric(d2);   qDiff <- numeric(d2);   qNew <- numeric(d2)
        p <- numeric(d2);   pL <- numeric(d2);   pR <- numeric(d2);   p2 <- numeric(d2);      p3 <- numeric(d2)
        grad <- numeric(d2);   gradFirst <- numeric(d2);   gradSaveL <- numeric(d2);   gradSaveR <- numeric(d2)
        log2 <- log(2)
        warningCodes <- array(0, c(max(numWarnings,1), 2))
        warningInd <- 0
        epsilonOrig <- epsilon
        warmupIntervalLengths <- rep(0,2)
        warmupIntervalsAdaptM <- rep(0,2)
        warmupIntervalNumber <- 0
        warmupIntervalCount <- 0
        epsilonAdaptCount <- 0
        warmupSamples <- array(0, c(2,d2))           ## 2xd array
        if(length(M) == 1) { if(M == -1) M <- rep(1, d2) else M <- c(M, 1) }
        Morig <- M
        sqrtM <- sqrt(M)
        numDivergences <- 0
        numTimesMaxTreeDepth <- 0
        ## nimbleLists
        qpNLDef <- nimbleList(q  = double(1), p  = double(1))
        btNLDef <- nimbleList(q1 = double(1), p1 = double(1), q2 = double(1), p2 = double(1), q3 = double(1), n = double(), s = double(), a = double(), na = double())
        ## nested function and function list definitions
        my_setWarmup <- hmc_setWarmup(warmupMode, warmup, messages, 'NUTS_classic', targetNodesToPrint)
        ## checks
        if(!isTRUE(nimbleOptions('enableDerivs')))   stop('must enable NIMBLE derivatives, set nimbleOptions(enableDerivs = TRUE)', call. = FALSE)
        if(!isTRUE(model$modelDef[['buildDerivs']])) stop('must set buildDerivs = TRUE when building model',  call. = FALSE)
        if(epsilon < 0) stop('NUTS_classic sampler epsilon must be non-negative', call. = FALSE)
        if(!all(M > 0)) stop('NUTS_classic sampler M must contain all positive elements', call. = FALSE)
        if(d == 1) if(length(M) != 2) stop('length of NUTS_classic sampler M must match length of NUTS_classic target nodes', call. = FALSE)
        if(d  > 1) if(length(M) != d) stop('length of NUTS_classic sampler M must match length of NUTS_classic target nodes', call. = FALSE)
        if(maxTreeDepth < 1) stop('NUTS_classic maxTreeDepth must be at least one ', call. = FALSE)
        hmc_checkWarmup(warmupMode, warmup, 'NUTS_classic')
    },
    run = function() {
        ## No-U-Turn Sampler with Dual Averaging, Algorithm 6 from Hoffman and Gelman (2014)
        if(timesRan == 0) {
            if(nwarmup == -1) stop('NUTS_classic nwarmup was not set correctly')
            ## reduce all pre-allocated vectors to correct size (d)
            q <<- q[1:d];   qL <<- qL[1:d];   qR <<- qR[1:d];   qDiff <<- qDiff[1:d];   qNew <<- qNew[1:d]
            p <<- p[1:d];   pL <<- pL[1:d];   pR <<- pR[1:d];   p2 <<- p2[1:d];           p3 <<- p3[1:d]
            grad <<- grad[1:d];   gradFirst <<- gradFirst[1:d];   gradSaveL <<- gradSaveL[1:d];   gradSaveR <<- gradSaveR[1:d]
            M <<- M[1:d];         sqrtM <<- sqrtM[1:d]
            if(epsilon == 0) epsilon <<- 1
            mu <<- log(10*epsilon)              ## following Stan: use default 1 and set mu before initializeEpsilon for first window
            if(initializeEpsilon & adaptive)  initEpsilon()
        }
        timesRan <<- timesRan + 1
        if(printTimesRan) print('============ times ran = ', timesRan)
        if(printEpsilon)  print('epsilon = ', epsilon)
        if(printM)        { print('M:'); print(M) }
        q <<- my_parameterTransform$transform(values(model, targetNodes))
        drawMomentumValues()    ## draws values for p
        qpLogH <- logH(q, p)
        logu <- qpLogH - rexp(1, 1)    ## logu <- lp - rexp(1, 1) => exp(logu) ~ uniform(0, exp(lp))
        qL <<- q;   qR <<- q;   pL <<- p;   pR <<- p;   j  <- 1;   n <- 1;   s <- 1;   qNew <<- q
        while(s == 1) {
            v <- 2*rbinom(1, 1, 0.5) - 1    ## -1 or 1
            if(v == -1) {
                btNL <- buildtree(qL, pL, logu, v, j, epsilon, qpLogH, 1)   ## first call: first = 1
                qL <<- btNL$q1;   pL <<- btNL$p1
            } else {
                btNL <- buildtree(qR, pR, logu, v, j, epsilon, qpLogH, 1)   ## first call: first = 1
                qR <<- btNL$q2;   pR <<- btNL$p2
            }
            if(btNL$s == 1)   if(runif(1) < btNL$n / n)   qNew <<- btNL$q3
            n <- n + btNL$n
            qDiff <<- qR - qL
            ##s <- btNL$s * nimStep(inprod(qDiff, pL)) * nimStep(inprod(qDiff, pR))                      ## this line replaced with the next,
            if(btNL$s == 0) s <- 0 else s <- nimStep(inprod(qDiff, pL)) * nimStep(inprod(qDiff, pR))     ## which acccounts for NaN's in btNL elements
            if(j >= maxTreeDepth) s <- 0
            if(printJ) {
                if(j == 1) cat('j = ', j) else cat(', ', j)
                cat('(');   if(v==1) cat('R') else cat('L');   cat(')')
                if(s != 1) print(' ')
            }
            if(j >= maxTreeDepth) { numTimesMaxTreeDepth <<- numTimesMaxTreeDepth + 1 }
            j <- j + 1
            checkInterrupt()
        }
        inverseTransformStoreCalculate(qNew)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        if((timesRan <= nwarmup) & adaptive)   adaptiveProcedure(btNL$a, btNL$na)
    },
    methods = list(
        drawMomentumValues = function() {
            ## M holds the diagonal elements of the diagonal mass matrix (aka the metric), which is:
            ## - the covariance matrix used for drawing p,
            ## - used in the calculation of kinetic energy K(p), and
            ## - used in the middle step of the leapfrog update.
            ## We let M approximate the precision matrix of q.
            for(i in 1:d)   p[i] <<- rnorm(1, 0, sqrtM[i])
        },
        logH = function(qArg = double(1), pArg = double(1)) {
            ## see comments in drawMomentumValues method
            lp <- inverseTransformStoreCalculate(qArg) - sum(pArg^2/M)/2 + my_parameterTransform$logDetJacobian(qArg)
            returnType(double())
            return(lp)
        },
        inverseTransformStoreCalculate = function(qArg = double(1)) {
            values(model, targetNodes) <<- my_parameterTransform$inverseTransform(qArg)
            lp <- model$calculate(calcNodes)
            returnType(double())
            return(lp)
        },
        calcLogProb = function(qArg = double(1)) {
            ans <- inverseTransformStoreCalculate(qArg) + my_parameterTransform$logDetJacobian(qArg)     ## won't forget this
            returnType(double())
            return(ans)
        },
        gradient_aux = function(qArg = double(1)) {
            derivsOutput <- nimDerivs(calcLogProb(qArg), order = 1, wrt = nimDerivs_wrt, model = model, updateNodes = nimDerivs_updateNodes, constantNodes = nimDerivs_constantNodes)
            returnType(double(1))
            return(derivsOutput$jacobian[1, 1:d])
        },
        gradient = function(qArg = double(1)) {
            derivsOutput <- nimDerivs(gradient_aux(qArg), order = 0, wrt = nimDerivs_wrt, model = model, updateNodes = nimDerivs_updateNodes, constantNodes = nimDerivs_constantNodes)
            grad <<- derivsOutput$value
        },
        leapfrog = function(qArg = double(1), pArg = double(1), eps = double(), first = double(), v = double()) {
            ## Algorithm 1 from Hoffman and Gelman (2014)
            if(first == 1) {
                gradient(qArg)                  ## member data 'grad' is set in gradient() method
            } else {
                if(v ==  1) grad <<- gradSaveR
                if(v == -1) grad <<- gradSaveL
                if(v ==  2) grad <<- gradSaveL
            }
            p2 <<- pArg + eps/2 * grad
            q2 <-  qArg + eps   * p2/M          ## see comments in drawMomentumValues method
            gradFirst <<- grad
            gradient(q2)                        ## member data 'grad' is set in gradient() method
            p3 <<- p2   + eps/2 * grad
            if(first == 1) {
                if(v ==  1) { gradSaveL <<- gradFirst;   gradSaveR <<- grad }
                if(v == -1) { gradSaveR <<- gradFirst;   gradSaveL <<- grad }
                if(v ==  2) { gradSaveL <<- gradFirst                       }
            } else {
                if(v ==  1) gradSaveR <<- grad
                if(v == -1) gradSaveL <<- grad
            }
            if(warningInd < numWarnings) if(is.nan.vec(c(q2, p3))) { warningInd <<- warningInd + 1; warningCodes[warningInd,1] <<- 1; warningCodes[warningInd,2] <<- timesRan } ## message code 1: print('  [Warning] NUTS_classic sampler (nodes: ', targetNodesToPrint, ') encountered a NaN value in leapfrog routine, with timesRan = ', timesRan)
            returnType(qpNLDef());   return(qpNLDef$new(q = q2, p = p3))
        },
        initEpsilon = function() {
            ## Algorithm 4 from Hoffman and Gelman (2014)
            savedCalcNodeValues <- values(model, calcNodes)
            q <<- my_parameterTransform$transform(values(model, targetNodes))
            p <<- numeric(d)        ## keep, sets 'p' to size d on first iteration
            drawMomentumValues()    ## draws values for p
            qpNL <- leapfrog(q, p, epsilon, 1, 2)            ## v = 2 is a special case for initializeEpsilon routine
            while(is.nan.vec(qpNL$q) | is.nan.vec(qpNL$p)) {              ## my addition
                ##if(numWarnings > 0) { print('  [Warning] NUTS_classic sampler (nodes: ', targetNodesToPrint, ') encountered NaN while initializing step-size; recommend better initial values')
                ##                      print('            reducing initial step-size'); numWarnings <<- numWarnings - 1 }
                epsilon <<- epsilon / 2                                   ## my addition
                qpNL <- leapfrog(q, p, epsilon, 0, 2)                     ## my addition
            }                                                             ## my addition
            qpLogH <- logH(q, p)
            a <- 2*nimStep(exp(logH(qpNL$q, qpNL$p) - qpLogH) - 0.5) - 1
            if(warningInd < numWarnings) if(is.nan(a)) { warningInd <<- warningInd + 1; warningCodes[warningInd,1] <<- 2; warningCodes[warningInd,2] <<- timesRan } ## message code 2: print('  [Warning] NUTS_classic sampler (nodes: ', targetNodesToPrint, ') caught acceptance prob = NaN in initializeEpsilon routine')
            ## while(a * (logH(qpNL$q, qpNL$p) - qpLogH) > -a * log2) {   ## replaced by simplified expression:
            while(a * (logH(qpNL$q, qpNL$p) - qpLogH + log2) > 0) {
                epsilon <<- epsilon * 2^a
                qpNL <- leapfrog(q, p, epsilon, 0, 2)        ## v = 2 is a special case for initializeEpsilon routine
            }
            values(model, calcNodes) <<- savedCalcNodeValues
        },
        adaptiveProcedure = function(a = double(), na = double()) {
            ## adapt epsilon:
            ## this is the "Dual Averaging" part of Algorithm 6 from Hoffman and Gelman (2014)
            if(adaptEpsilon) {
                epsilonAdaptCount <<- epsilonAdaptCount + 1
                Hbar <<- (1 - 1/(epsilonAdaptCount+t0)) * Hbar + 1/(epsilonAdaptCount+t0) * (delta - a/na)
                logEpsilon <- mu - sqrt(epsilonAdaptCount)/gamma * Hbar
                epsilon <<- exp(logEpsilon)
                timesRanToNegativeKappa <- epsilonAdaptCount^(-kappa)
                logEpsilonBar <<- timesRanToNegativeKappa * logEpsilon + (1 - timesRanToNegativeKappa) * logEpsilonBar
                if(timesRan == nwarmup)   epsilon <<- exp(logEpsilonBar)
                if(warningInd < numWarnings) if(is.nan(epsilon)) { warningInd <<- warningInd + 1; warningCodes[warningInd,1] <<- 3; warningCodes[warningInd,2] <<- timesRan } ## message code 3: print('  [Warning] NUTS_classic sampler (nodes: ', targetNodesToPrint, ') value of epsilon is NaN, with timesRan = ', timesRan)
            }
            ## adapt M:
            if(adaptM) {
                if(warmupIntervalNumber <= length(warmupIntervalLengths)) {
                    warmupIntervalCount <<- warmupIntervalCount + 1
                    if(warmupIntervalCount > warmupIntervalLengths[warmupIntervalNumber]) stop('Unexpected behavior in NUTS_classic warmup book-keeping')
                    if(warmupIntervalsAdaptM[warmupIntervalNumber] == 1)   warmupSamples[warmupIntervalCount, 1:d] <<- qNew
                    if(warmupIntervalCount == warmupIntervalLengths[warmupIntervalNumber]) {
                        if(warmupIntervalsAdaptM[warmupIntervalNumber] == 1) {
                            ## see comments in drawMomentumValues method
                            ## use regularized estimation of empirical covariance (identical to Stan):
                            ## https://github.com/stan-dev/stan/blob/develop/src/stan/mcmc/covar_adaptation.hpp
                            ## SigmaRegularized = [ N / (N + 5) ] * Sigma_raw + 0.001 * [ 5 / (N + 5) ] * I
                            ## only estimating diagonal elements:
                            for(i in 1:d) {
                                v <- var(warmupSamples[1:warmupIntervalCount, i])
                                vReg <- (warmupIntervalCount/(warmupIntervalCount+5))*v + 0.001*(5/(warmupIntervalCount+5))
                                M[i] <<- 1/vReg
                            }
                            ## estimating full empirical covariance:
                            ##for(i in 1:d)     warmupSamples[, i] <- warmupSamples[, i] - mean(warmupSamples[, i])
                            ##warmupSamplesCov <- (t(warmupSamples) %*% warmupSamples) / (warmupIntervalCount-1)
                            ##warmupCovRegularized <- (warmupIntervalCount/(warmupIntervalCount+5))*warmupSamplesCov + 0.001*(5/(warmupIntervalCount+5))*diag(d)
                            ##for(i in 1:d)   M[i] <<- 1 / warmupCovRegularized[i,i]
                            sqrtM <<- sqrt(M)
                            if(adaptEpsilon) {
                                initEpsilon()
                                epsilonAdaptCount <<- 0
                                mu <<- log(10 * epsilon)
                            }
                        }
                        warmupIntervalCount <<- 0
                        warmupIntervalNumber <<- warmupIntervalNumber + 1
                    }
                }
            }
        },
        buildtree = function(qArg = double(1), pArg = double(1), logu = double(), v = double(), j = double(), eps = double(), logH0 = double(), first = double()) {
            ## Algorithm 6 (second half) from Hoffman and Gelman (2014)
            returnType(btNLDef())
            if(j == 1) {    ## one leapfrog step in the direction of v
                qpNL <- leapfrog(qArg, pArg, v*eps, first, v)
                q <<- qpNL$q;   p <<- qpNL$p;   qpLogH <- logH(q, p)
                n <- nimStep(qpLogH - logu)          ## step(x) = 1 iff x >= 0, and zero otherwise
                s <- nimStep(qpLogH - logu + deltaMax)
                ## lowering the initial step size, and increasing the target acceptance rate may keep the step size small to avoid divergent paths.
                if(s == 0) { numDivergences <<- numDivergences + 1 }
                ##           if(numWarnings > 0) { print('  [Warning] NUTS_classic sampler (nodes: ', targetNodesToPrint, ') encountered a divergent path on iteration ', timesRan, ', with divergence = ', logu - qpLogH)
                ##                                 numWarnings <<- numWarnings - 1 } }
                a <- min(1, exp(qpLogH - logH0))
                if(is.nan.vec(q) | is.nan.vec(p) | is.nan(a)) { n <- 0; s <- 0; a <- 0 }     ## my addition
                return(btNLDef$new(q1 = q, p1 = p, q2 = q, p2 = p, q3 = q, n = n, s = s, a = a, na = 1))
            } else {        ## recursively build left and right subtrees
                btNL1 <- buildtree(qArg, pArg, logu, v, j-1, eps, logH0, 0)
                if(btNL1$s == 1) {
                    if(v == -1) {
                        btNL2 <- buildtree(btNL1$q1, btNL1$p1, logu, v, j-1, eps, logH0, 0)   ## recursive calls: first = 0
                        btNL1$q1 <- btNL2$q1;   btNL1$p1 <- btNL2$p1
                    } else {
                        btNL2 <- buildtree(btNL1$q2, btNL1$p2, logu, v, j-1, eps, logH0, 0)   ## recursive calls: first = 0
                        btNL1$q2 <- btNL2$q2;   btNL1$p2 <- btNL2$p2
                    }
                    nSum <- btNL1$n + btNL2$n
                    if(nSum > 0)   if(runif(1) < btNL2$n / nSum)   btNL1$q3 <- btNL2$q3
                    qDiff <<- btNL1$q2 - btNL1$q1
                    btNL1$a  <- btNL1$a  + btNL2$a
                    btNL1$na <- btNL1$na + btNL2$na
                    btNL1$s  <- btNL2$s * nimStep(inprod(qDiff, btNL1$p1)) * nimStep(inprod(qDiff, btNL1$p2))
                    btNL1$n  <- nSum
                }
                return(btNL1)
            }
        },
        before_chain = function(MCMCniter = double(), MCMCnburnin = double(), MCMCchain = double()) {
            if(MCMCchain == 1)   nwarmup <<- my_setWarmup$run(MCMCniter, MCMCnburnin, adaptive)
            if(adaptive) {
                if(nwarmup > 0) {
                    ## need to deal with exceptions such as nwarmup < 100
                    if(initBuffer + adaptWindow + termBuffer > nwarmup) {
                        if(messages & adaptive) print('  [Warning] Number of warmup iterations for NUTS_classic sampler ',
                                                      'is too small for even one cycle of standard adaptation. Using 15% ',
                                                      'for initial stepsize adaptation, 75% for mass matrix and stepsize ',
                                                      'adaptation, and 10% for final stepsize adaptation.')
                        initBuffer <<- round(nwarmup * 0.15)
                        termBuffer <<- round(nwarmup * 0.10)
                        adaptWindow <<- nwarmup - initBuffer - termBuffer
                    }

                    warmupIntervalLengths <<- numeric(length = 1, value = initBuffer)
                    endIntervals <- initBuffer     ## iteration marking the end of intervals planned so far
                    warmupIntervalsAdaptM <<- numeric(length = 1, value = 0)
                    nextIntervalLength <- adaptWindow
                    done <- FALSE
                    while(!done) {
                        warmupIntervalLengths <<- c(warmupIntervalLengths, nextIntervalLength)
                        warmupIntervalsAdaptM <<- c(warmupIntervalsAdaptM, 1)
                        endIntervals <- endIntervals + nextIntervalLength
                        remainingIterations <- nwarmup - (endIntervals + termBuffer)
                        if(remainingIterations == 0) {
                            done <- TRUE
                        } else {
                            ## look ahead two iterations into the future
                            nextIntervalLength <- 2*nextIntervalLength
                            nextRemainingIterations <- remainingIterations - nextIntervalLength
                            nextnextIntervalLength <- 2*nextIntervalLength
                            if(nextRemainingIterations < nextnextIntervalLength) {
                                nextIntervalLength <- nextIntervalLength + nextRemainingIterations
                            }
                        }
                    }
                    warmupIntervalLengths <<- c(warmupIntervalLengths, termBuffer)
                    warmupIntervalsAdaptM <<- c(warmupIntervalsAdaptM, 0)
                    setSize(warmupSamples, max(warmupIntervalLengths), d, fillZeros = FALSE)
                }
            }
        },
        after_chain = function() {
            if(messages) {
                if(numDivergences == 1) print('  [Note] NUTS_classic sampler (nodes: ', targetNodesToPrint, ') encountered ', numDivergences, ' divergent path.')
                if(numDivergences  > 1) print('  [Note] NUTS_classic sampler (nodes: ', targetNodesToPrint, ') encountered ', numDivergences, ' divergent paths.')
                if(numTimesMaxTreeDepth == 1) print('  [Note] NUTS_classic sampler (nodes: ', targetNodesToPrint, ') reached the maximum search tree depth ', numTimesMaxTreeDepth, ' time.')
                if(numTimesMaxTreeDepth  > 1) print('  [Note] NUTS_classic sampler (nodes: ', targetNodesToPrint, ') reached the maximum search tree depth ', numTimesMaxTreeDepth, ' times.')
                numDivergences <<- 0           ## reset counters for numDivergences and numTimesMaxTreeDepth,
                numTimesMaxTreeDepth <<- 0     ## even when using reset=FALSE to continue the same chain
            }
            if(warningInd > 0) {
                for(i in 1:warningInd) {
                    if(warningCodes[i,1] == 1) print('  [Warning] NUTS_classic sampler (nodes: ', targetNodesToPrint, ') encountered a NaN value on MCMC iteration ', warningCodes[i,2], '.')
                    if(warningCodes[i,1] == 2) print('  [Warning] NUTS_classic sampler (nodes: ', targetNodesToPrint, ') encountered acceptance prob = NaN in initEpsilon routine.')
                    if(warningCodes[i,1] == 3) print('  [Warning] NUTS_classic sampler (nodes: ', targetNodesToPrint, ') encountered epsilon = NaN on MCMC iteration ', warningCodes[i,2], '.')
                }
                warningInd <<- 0               ## reset warningInd even when using reset=FALSE to continue the same chain
            }
        },
        reset = function() {
            timesRan       <<- 0
            epsilon        <<- epsilonOrig
            mu             <<- 0
            logEpsilonBar  <<- 0
            Hbar           <<- 0
            numDivergences <<- 0
            numTimesMaxTreeDepth <<- 0
            warningInd     <<- 0
            M              <<- Morig
            sqrtM          <<- sqrt(M)
            warmupIntervalNumber <<- 1
            warmupIntervalCount  <<- 0
        }
    ),
    buildDerivs = list(
        inverseTransformStoreCalculate = list(),
        calcLogProb = list(),
        gradient_aux = list()
    )
)



#' nimbleList definition used internally in NUTS sampler.
#' @export
stateNL_NUTS <- nimbleList(q = double(1), p = double(1), H = double(), logProb = double(), gr_logProb = double(1))

#' nimbleList definition used internally in NUTS sampler.
#' @export
treebranchNL_NUTS <- nimbleList(p_beg = double(1), p_end = double(1), rho = double(1), log_sum_wt = double())

#' No-U-Turn (NUTS) Hamiltonian Monte Carlo (HMC) Sampler
#'
#' The NUTS sampler implements No-U-Turn (NUTS) Hamiltonian Monte Carlo (HMC) sampling following the algorithm of version 2.32.2 of Stan. Internally, any posterior dimensions with bounded support are transformed, so sampling takes place on an unconstrained space. In contrast to standard HMC (Neal, 2011), the NUTS algorithm removes the tuning parameters of the leapfrog step size and the number of leapfrog steps, thus providing a sampling algorithm that can be used without hand-tuning or trial runs.
#'
#' @param model An uncompiled nimble model object on which the MCMC will operate.
#' @param mvSaved A nimble \code{modelValues} object to be used to store MCMC samples.
#' @param target A character vector of node names on which the sampler will operate.
#' @param control A named list that controls the precise behavior of the sampler. The default values for control list elements are specified in the setup code of the sampler. A description of the possible control list elements appear in the details section.
#'
#' @details
#'
#' The NUTS sampler accepts the following control list elements:
#' 
#' \itemize{
#' \item messages. A logical argument, specifying whether to print informative messages (default = TRUE)
#' \item numWarnings. A numeric argument, specifying how many warnings messages to emit (for example, when \code{NaN} values are encountered). See additional details below. (default = 0)
#' \item epsilon. A positive numeric argument, specifying the initial step-size value. If not provided, an appropriate initial value is selected.
#' \item gamma. A positive numeric argument, specifying the degree of shrinkage used during the initial period of step-size adaptation. (default = 0.05)
#' \item t0. A non-negative numeric argument, where larger values stabilize (attenuate) the initial period of step-size adaptation. (default = 10)
#' \item kappa. A numeric argument between zero and one, where smaller values give a higher weighting to more recent iterations during the initial period of step-size adaptation. (default = 0.75)
#' \item delta. A numeric argument, specifying the target acceptance probability used during the initial period of step-size adaptation. (default = 0.8)
#' \item deltaMax. A positive numeric argument, specifying the maximum allowable divergence from the Hamiltonian value. Paths which exceed this value are considered divergent, and will not proceed further. (default = 1000)
#' \item M. A vector of positive real numbers, with length equal to the number of dimensions being sampled. Elements of \code{M} specify the diagonal elements of the diagonal mass matrix (or the metric) used for the auxiliary momentum variables in sampling. Sampling may be improved if the elements of \code{M} approximate the marginal inverse variance (precision) of the (potentially transformed) parameters. (default: a vector of ones).
#' \item warmupMode. A character string, specifying the behavior for choosing the number of warmup iterations. Four values are possible. The value 'default' (the default) sets the number of warmup iterations as the number of burnin iterations (if a positive value for \code{nburnin} is used) or half the number of MCMC iterations in each chain (if \code{nburnin = 0}). The value 'burnin' sets the number of warmup iterations as the number of burnin iterations regardless of the length of the burnin period. The value 'fraction' sets the number of warmup iterations as \code{fraction*niter}, where \code{fraction} is the value of the \code{warmup} control argument, and \code{niter} is the number of MCMC iterations in each chain; in this case, the value of the \code{warmup} control argument must be between 0 and 1. The value 'iterations' sets the number of warmup iterations as the value of the \code{warmup} control argumnet, regardless of the length of the burnin period or the number of MCMC iterations; in this case the value of \code{warmup} must be a non-negative integer. In all cases, the number of (pre-thinning) samples discarded equals \code{nburnin}, as is always the case for MCMC in NIMBLE.
#' \item warmup. Numeric value used in determining the number of warmup iterations. This control argument is only used when \code{warmupMode} is 'fraction' or 'iterations'. 
#' \item maxTreeDepth. The maximum allowable depth of the binary leapfrog search tree for generating candidate transitions. (default = 10)
#' \item adaptWindow. Number of iterations in the first adaptation window used for adapting the mass matrix (M). Subsequent adaptation windows double in length, so long as enough warmup iterations are available. (default = 25)
#' \item initBuffer. Number of iterations in the initial warmup window, which occurs prior to the first adaptation of the metric M. (default = 75)
#' \item termBuffer. Number of iterations in the final (terminal) warmup window, before which the metric M is not adjusted(default = 50)
#' \item adaptive. A logical argument, specifying whether to do any adaptation whatsoever. When \code{TRUE}, specific adaptation routines are controlled by the \code{adaptEpsilon} and \code{adaptM} control list elements. (default = TRUE)
#' \item adaptEpsilon. A logical argument, specifying whether to perform stepsize adaptation. Only used when \code{adaptive = TRUE}. (default = TRUE)
#' \item adaptM. A logical argument, specifying whether to perform adaptation of the mass matrix (metric) M. Only used when \code{adaptive = TRUE}. (default = TRUE)
#' \item initializeEpsilon. A logical argument, specifying whether to perform the epsilon (stepsize) initialization routine at the onset of each adaptation window. (default = TRUE)
#' }
#'
#' \code{NaN} values may be encountered in the course of the leapfrog procedure. In particular, when the stepsize (epsilon) is too large, the leapfrog procedure can step too far and arrive at an invalid region of parameter space, thus generating a \code{NaN} value in the likelihood evaluation or in the gradient calculation. These situation are handled by the sampler by rejecting the \code{NaN} value, and reducing the stepsize.
#' 
#' @import nimble
#' 
#' @export
#'
#' @return A object of class `sampler_NUTS`.
#' 
#' @aliases NUTS nuts HMC hmc sampler_NUTS
#' 
#' @author Perry de Valpine and Daniel Turek
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
#' conf <- configureMCMC(Rmodel, nodes = NULL)
#' 
#' conf$addSampler(target = c('b0', 'b1', 'sigma'), type = 'NUTS')
#' 
#' Rmcmc <- buildMCMC(conf)
#'
#' @references
#'
#' Hoffman, Matthew D., and Gelman, Andrew (2014). The No-U-Turn Sampler: Adaptively setting path lengths in Hamiltonian Monte Carlo. \emph{Journal of Machine Learning Research}, 15(1): 1593-1623.
#'
#' Stan Development Team. 2023. Stan Modeling Language Users Guide and Reference Manual, 2.32.2. https://mc-stan.org.
sampler_NUTS <- nimbleFunction(
    name = 'sampler_NUTS',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        printTimesRan  <- extractControlElement(control, 'printTimesRan',  FALSE)
        printEpsilon   <- extractControlElement(control, 'printEpsilon',   FALSE)
        printJ         <- extractControlElement(control, 'printJ',         FALSE)
        printM         <- extractControlElement(control, 'printM',         FALSE)
        messages       <- extractControlElement(control, 'messages',       getNimbleOption('verbose'))
        numWarnings    <- extractControlElement(control, 'numWarnings',    0)
        epsilon        <- extractControlElement(control, 'epsilon',        0)          # initial epsilon, if 0 then use 1
        initializeEpsilon <- extractControlElement(control, 'initializeEpsilon', TRUE) # use initializeEpsilon step?
        gamma          <- extractControlElement(control, 'gamma',          0.05)
        t0             <- extractControlElement(control, 't0',             10)
        kappa          <- extractControlElement(control, 'kappa',          0.75)
        delta          <- extractControlElement(control, 'delta',          0.8)
        deltaMax       <- extractControlElement(control, 'deltaMax',       1000)
        M              <- extractControlElement(control, 'M',              -1)
        warmupMode     <- extractControlElement(control, 'warmupMode',     'default')   ## 'default', 'burnin', 'fraction', or 'iterations'
        warmup         <- extractControlElement(control, 'warmup',         -1)          ## used if warmupMode is 'fraction' or 'iterations'
        maxTreeDepth   <- extractControlElement(control, 'maxTreeDepth',   10)
        adaptWindow    <- extractControlElement(control, 'adaptWindow',    25)
        initBuffer     <- extractControlElement(control, 'initBuffer',     75)
        termBuffer     <- extractControlElement(control, 'termBuffer',     50)
        adaptive       <- extractControlElement(control, 'adaptive',       TRUE) # any adaptation? (if FALSE, next two flags are ignored)
        adaptEpsilon   <- extractControlElement(control, 'adaptEpsilon',   TRUE) # stepsize adaptation?
        adaptM         <- extractControlElement(control, 'adaptM',         TRUE) # mass matrix adaptation?
        ## node list generation
        targetNodes <- model$expandNodeNames(target)
        if(length(targetNodes) <= 0) stop('NUTS sampler must operate on at least one node', call. = FALSE)
        targetNodesAsScalars <- model$expandNodeNames(targetNodes, returnScalarComponents = TRUE)
        targetNodesToPrint <- paste(targetNodes, collapse = ', ')
        if(nchar(targetNodesToPrint) > 100)   targetNodesToPrint <- paste0(substr(targetNodesToPrint, 1, 97), '...')
        calcNodes <- model$getDependencies(targetNodes)
        ## check for discrete nodes (early, before parameterTransform is specialized)
        if(any(model$isDiscrete(targetNodesAsScalars)))
            stop(paste0('NUTS sampler cannot operate on discrete-valued nodes: ',
                        paste0(targetNodesAsScalars[model$isDiscrete(targetNodesAsScalars)], collapse = ', ')))
        ## processing of bounds and transformations
        my_parameterTransform <- parameterTransform(model, targetNodesAsScalars)
        d <- my_parameterTransform$getTransformedLength()
        d2 <- max(d, 2) ## for pre-allocating vectors
        nimDerivs_wrt <- 1:d
        derivsInfo_return <- makeModelDerivsInfo(model, targetNodes, calcNodes)
        nimDerivs_updateNodes   <- derivsInfo_return$updateNodes
        nimDerivs_constantNodes <- derivsInfo_return$constantNodes
        ## numeric value generation
        epsilonOrig <- epsilon
        timesRan            <- 0
        nwarmup             <- 0
        stepsizeCounter     <- 0
        mu                  <- 0
        logEpsilonBar       <- 0
        Hbar                <- 0
        n_leapfrog          <- 0
        sum_metropolis_prob <- 0
        warmupSamples       <- array(0, c(2,d2))     ## 2xd array
        divergent           <- FALSE
        log2                <- log(2)
        ## adapt_* variables are initialized in before_chain method
        adaptWindow_size    <- 0
        adapt_initBuffer    <- 0
        adapt_termBuffer    <- 0
        adapt_next_window   <- 0
        adaptWindow_counter <- 0
        adaptWindow_iter    <- 0
        if(length(M) == 1) { if(M == -1) M <- rep(1, d2) else M <- c(M, 1) }
        Morig <- M
        sqrtM <- sqrt(M)
        warningInd   <- 0
        warningCodes <- array(0, c(max(numWarnings,1), 2))
        numDivergences       <- 0
        numTimesMaxTreeDepth <- 0
        ## nimbleLists
        treebranchNL <- treebranchNL_NUTS   ## reference input to buildtree
        stateNL <- stateNL_NUTS             ## system state (p, q, H, lp, gr_lp)
        state_current <- stateNL$new()
        state_f       <- stateNL$new()
        state_b       <- stateNL$new()
        state_sample  <- stateNL$new()
        state_propose <- stateNL$new()
        ## nested function and function list definitions
        my_setWarmup <- hmc_setWarmup(warmupMode, warmup, messages, 'NUTS', targetNodesToPrint)
        ## checks
        if(!isTRUE(nimbleOptions('enableDerivs')))   stop('must enable NIMBLE derivatives, set nimbleOptions(enableDerivs = TRUE)', call. = FALSE)
        if(!isTRUE(model$modelDef[['buildDerivs']])) stop('must set buildDerivs = TRUE when building model',  call. = FALSE)
        if(epsilon < 0) stop('NUTS sampler epsilon must be non-negative', call. = FALSE)
        if(!all(M > 0)) stop('NUTS sampler M must contain all positive elements', call. = FALSE)
        if(d == 1) if(length(M) != 2) stop('length of NUTS sampler M must match length of NUTS target nodes', call. = FALSE)
        if(d  > 1) if(length(M) != d) stop('length of NUTS sampler M must match length of NUTS target nodes', call. = FALSE)
        if(maxTreeDepth < 1) stop('NUTS maxTreeDepth must be at least one', call. = FALSE)
        hmc_checkWarmup(warmupMode, warmup, 'NUTS')
    },
    run = function() {
        ## No-U-Turn Sampler based on Stan
        state_current$q <<- my_parameterTransform$transform(values(model, targetNodes))
        if(timesRan == 0) {
            if(nwarmup == -1) stop('NUTS nwarmup was not set correctly')
            state_current$p          <<- numeric(d, init = FALSE)
            state_current$gr_logProb <<- numeric(d, init = FALSE)
            M <<- M[1:d]
            sqrtM <<- sqrtM[1:d]
            if(epsilon <= 0) epsilon <<- 1
            mu <<- log(10*epsilon)    ## curiously, Stan sets this for the first round *before* init_stepsize
            if(initializeEpsilon & adaptive)   initEpsilon()
        }
        timesRan <<- timesRan + 1
        if(printTimesRan) print('============ times ran = ', timesRan)
        if(printEpsilon)  print('epsilon = ', epsilon)
        if(printM)        { print('M:'); print(M) }
        drawMomentumValues(state_current)    ## draws values for p
        update_state_calcs(state_current)
        copy_state(state_f, state_current)
        copy_state(state_b, state_current)
        copy_state(state_sample, state_current)
        copy_state(state_propose, state_current)
        ##
        p_ff <- state_current$p
        p_fb <- state_current$p
        p_bf <- state_current$p
        p_bb <- state_current$p
        ##
        rho <- state_current$p
        ##
        log_sum_wt <- 0
        H0 <- state_current$H
        depth <- 0
        n_leapfrog <<- 0
        sum_metropolis_prob <<- 0
        divergent <<- FALSE
        branch <- treebranchNL$new()
        done <- FALSE
        old_stepsize <- epsilon
        old_M <- M
        while((depth < maxTreeDepth) & (!done)) {
            checkInterrupt()
            rho_f <- numeric(0, length = d)
            rho_b <- numeric(0, length = d)
            valid_subtree <- FALSE
            log_sum_wt_subtree <- -Inf
            bool_fwd <- runif(1,0,1) > 0.5
            if(bool_fwd) {
                copy_state(state_current, state_f)
                rho_b <- rho
                p_bf <- p_ff
                ##
                branch$p_beg <- p_fb
                branch$p_end <- p_ff
                branch$rho <- rho_f
                branch$log_sum_wt <- log_sum_wt_subtree
                valid_subtree <- buildtree(depth, state_propose, branch, H0, 1)
                log_sum_wt_subtree <- branch$log_sum_wt
                rho_f <- branch$rho
                p_ff <- branch$p_end
                p_fb <- branch$p_beg
                copy_state(state_f, state_current)
            } else {
                copy_state(state_current, state_b)
                rho_f <- rho
                p_fb <- p_bb
                branch$p_beg <- p_bf
                branch$p_end <- p_bb
                branch$rho <- rho_b
                branch$log_sum_wt <- log_sum_wt_subtree
                valid_subtree <- buildtree(depth, state_propose, branch, H0, -1)
                log_sum_wt_subtree <- branch$log_sum_wt
                rho_b <- branch$rho
                p_bb <- branch$p_end
                p_bf <- branch$p_beg
                copy_state(state_b, state_current)
            }
            if(!valid_subtree)   done <- TRUE
            if(!done) {
                depth <- depth + 1
                accept <- FALSE
                log_accept_prob <- log_sum_wt_subtree - log_sum_wt
                if(log_accept_prob > 0) {
                    accept <- TRUE
                } else {
                    if(runif(1,0,1) < exp(log_accept_prob))   accept <- TRUE
                }
                if(accept)   copy_state(state_sample, state_propose)
                log_sum_wt <- log_sum_exp(log_sum_wt, log_sum_wt_subtree)
                rho <- rho_b + rho_f
                ##
                persist_criterion <- decide_persist(p_bb, p_ff, p_bf, p_fb, rho, rho_f, rho_b)
                if(!persist_criterion)   done <- TRUE
            }
        }
        ##
        accept_prob <- sum_metropolis_prob / n_leapfrog
        copy_state(state_current, state_sample)        ## extraneous copy? could remove?
        ##
        inverseTransformStoreCalculate(state_sample$q)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        if((timesRan <= nwarmup) & adaptive) {
            if(adaptEpsilon)   adapt_stepsize(accept_prob)
            update <- FALSE
            if(adaptM)   update <- adapt_M()
            if(update & adaptEpsilon) {
                if(initializeEpsilon)   initEpsilon()
                Hbar <<- 0
                logEpsilonBar <<- 0
                stepsizeCounter <<- 0
                mu <<- log(10*epsilon)
            }
        }
    },
    methods = list(
        copy_state = function(to = stateNL(), from = stateNL()) {
            to$p          <- from$p
            to$q          <- from$q
            to$H          <- from$H
            to$logProb    <- from$logProb
            to$gr_logProb <- from$gr_logProb
        },
        drawMomentumValues = function(state = stateNL()) {
            draws <- rnorm(length(sqrtM), 0, sd = sqrtM)
            for(i in 1:d)   state$p[i] <- draws[i]
        },
        inverseTransformStoreCalculate = function(qArg = double(1)) {
            values(model, targetNodes) <<- my_parameterTransform$inverseTransform(qArg)
            lp <- model$calculate(calcNodes)
            returnType(double())
            return(lp)
        },
        calcLogProb = function(qArg = double(1)) {
            ans <- inverseTransformStoreCalculate(qArg) + my_parameterTransform$logDetJacobian(qArg)     ## won't forget this
            returnType(double())
            return(ans)
        },
        gradient_aux = function(qArg = double(1)) {
            derivsOutput <- nimDerivs(calcLogProb(qArg), order = 1, wrt = nimDerivs_wrt, model = model, updateNodes = nimDerivs_updateNodes, constantNodes = nimDerivs_constantNodes)
            returnType(double(1))
            return(derivsOutput$jacobian[1, 1:d])
        },
        gradient = function(qArg = double(1)) {
            derivsOutput <- nimDerivs(gradient_aux(qArg), order = 0, wrt = nimDerivs_wrt, model = model, updateNodes = nimDerivs_updateNodes, constantNodes = nimDerivs_constantNodes)
            returnType(double(1))
            return(derivsOutput$value)
        },
        update_state_calcs = function(state = stateNL()) {
            state$logProb <- inverseTransformStoreCalculate(state$q) + my_parameterTransform$logDetJacobian(state$q)
            state$H <- 0.5 * sum(state$p * state$p/M) - state$logProb
            state$gr_logProb <- gradient(state$q)
        },
        leapfrog = function(state = stateNL(), eps = double()) {
            state$p <- state$p + 0.5 * eps * state$gr_logProb
            state$q <- state$q +       eps * state$p / M
            update_state_calcs(state)
            state$p <- state$p + 0.5 * eps * state$gr_logProb
            ## update H for the closing step in p that is after update_state_calcs
            state$H <- 0.5 * sum(state$p * state$p / M) - state$logProb
        },
        log_sum_exp = function(x1 = double(), x2 = double()) {
            C <- max(x1, x2)
            ans <- C + log(exp(x1 - C) + exp(x2 - C))
            returnType(double())
            return(ans)
        },
        buildtree = function(depth = integer(), state_propose = stateNL(), branch = treebranchNL(), H0 = double(), sign = double()) {
            if(depth == 0) {
                leapfrog(state_current, sign*epsilon)
                n_leapfrog <<- n_leapfrog + 1
                new_H <- state_current$H
                if(is.nan(new_H))   new_H <- Inf
                deltaH <- new_H - H0
                if(deltaH > deltaMax)   divergent <<- TRUE
                branch$log_sum_wt <- log_sum_exp(branch$log_sum_wt, -deltaH)
                if((-deltaH) > 0)   sum_metropolis_prob <<- sum_metropolis_prob + 1
                else                sum_metropolis_prob <<- sum_metropolis_prob + exp(-deltaH)
                copy_state(state_propose, state_current)
                branch$rho <- branch$rho + state_current$p
                branch$p_beg <- state_current$p
                branch$p_end <- state_current$p
                return(!divergent)
            }
            init_branch <- treebranchNL$new()
            init_branch$log_sum_wt <- -Inf
            init_branch$p_beg <- branch$p_beg
            init_branch$rho <- numeric(0, length = d)
            valid_init <- buildtree(depth-1, state_propose, init_branch, H0, sign)
            branch$p_beg <- init_branch$p_beg
            if(!valid_init)   return(FALSE)
            ##
            state_propose_final <- stateNL$new()
            copy_state(state_propose_final, state_current)
            final_branch <- treebranchNL$new()
            final_branch$p_end <- branch$p_end
            final_branch$log_sum_wt <- -Inf
            final_branch$rho <- numeric(0, length = d)
            valid_final <- buildtree(depth-1, state_propose_final, final_branch, H0, sign)
            branch$p_end <- final_branch$p_end
            if(!valid_final)   return(FALSE)
            ##
            log_sum_wt_subtree <- log_sum_exp(init_branch$log_sum_wt, final_branch$log_sum_wt)
            branch$log_sum_wt <- log_sum_exp(branch$log_sum_wt, log_sum_wt_subtree)
            ##
            accept <- FALSE
            log_accept_prob <- final_branch$log_sum_wt - log_sum_wt_subtree
            if(log_accept_prob > 0) {
                accept <- TRUE
            } else {
                accept_prob <- exp(log_accept_prob)
                accept <- runif(1, 0, 1) < accept_prob
            }
            if(accept)   copy_state(state_propose, state_propose_final)
            ##
            rho_subtree <- final_branch$rho + init_branch$rho
            branch$rho <- branch$rho + rho_subtree
            ##
            persist_criterion <- decide_persist(branch$p_beg, branch$p_end, init_branch$p_end,final_branch$p_beg, rho_subtree, final_branch$rho, init_branch$rho)
            ## For the case of diagonal M, the "sharp" variables in Stan are simply element-wise division by diag(M), e.g.:
            ## p_sharp_beg <- branch$p_beg / M
            ## p_sharp_end <- branch$p_end / M
            ## but these are now in the decide_persist function
            returnType(logical())
            return(persist_criterion)
        },
        decide_persist = function(p_b = double(1), p_e = double(1), p_b2 = double(1), p_e2 = double(1), rho = double(1), rho_1 = double(1), rho_2 = double(1)) {
            p_b_overM <- p_b / M
            p_e_overM <- p_e / M
            persist_criterion <- compute_criterion(p_b_overM, p_e_overM, rho)
            if(persist_criterion) {
                rho_alt <- rho_2 + p_e2
                p_e2_overM <- p_e2 / M
                persist_criterion <- compute_criterion(p_b_overM, p_e2_overM, rho_alt)
            }
            if(persist_criterion) {
                rho_alt <- rho_1 + p_b2
                p_b2_overM <- p_b2 / M
                persist_criterion <- compute_criterion(p_b2_overM, p_e_overM, rho_alt)
            }
            returnType(logical())
            return(persist_criterion)
        },
        compute_criterion = function(poverM1 = double(1), poverM2 = double(1), rho = double(1)) {
            ans <- (inprod(poverM2, rho) > 0) & (inprod(poverM1, rho) > 0)
            returnType(logical())
            return(ans)
        },
        initEpsilon = function() {
            initValues <- values(model, calcNodes)
            state_init <- stateNL$new()
            copy_state(state_init, state_current)
            drawMomentumValues(state_current)    ## draws values for p
            update_state_calcs(state_current)
            H0 <- state_current$H
            leapfrog(state_current, epsilon)
            newH <- state_current$H
            if(is.nan(newH))   newH <- Inf
            deltaH <- H0 - newH
            direction <- 2*nimStep(deltaH > log(0.8)) - 1    ## 1 or -1
            done <- FALSE
            while(!done) {
                copy_state(state_current, state_init)
                drawMomentumValues(state_current)
                update_state_calcs(state_current)
                H0 <- state_current$H
                leapfrog(state_current, epsilon)
                newH <- state_current$H
                if(is.nan(newH))   newH <- Inf
                deltaH <- H0 - newH
                if((direction == 1) & !(deltaH > log(0.8)))   done <- TRUE
                else {
                    if((direction == -1) & !(deltaH < log(0.8)))   done <- TRUE
                    else {
                        if(direction == 1) epsilon <<- epsilon * 2
                        else               epsilon <<- epsilon * 0.5
                    }
                }
                if(!done) {
                    if(epsilon > 1e7)    stop("Search for initial stepsize in NUTS sampler exploded. Something is wrong.")
                    if(epsilon < 1e-16)  stop("Search for initial stepsize in NUTS sampler shrank to 0. Something is wrong.")
                }
            }
            copy_state(state_current, state_init)
            values(model, calcNodes) <<- initValues
        },
        adapt_stepsize = function(adapt_stat = double()) {
            ## following Stan code, this is the same as what we have from Hoffman and Gelman, but with adapt_stat instead of a/na
            if(adapt_stat > 1)   adapt_stat <- 1
            old_epsilon <- epsilon
            stepsizeCounter <<- stepsizeCounter + 1
            eta <- 1/(stepsizeCounter + t0)
            Hbar <<- (1-eta) * Hbar + eta * (delta - adapt_stat)          ## s_bar in Stan
            logEpsilon <- mu - Hbar * sqrt(stepsizeCounter)/gamma         ## x in Stan
            epsilon <<- exp(logEpsilon)
            stepsizeCounterToNegativeKappa <- stepsizeCounter^(-kappa)    ## x_eta in Stan
            logEpsilonBar <<- stepsizeCounterToNegativeKappa * logEpsilon + (1 - stepsizeCounterToNegativeKappa) * logEpsilonBar   ## x_bar in Stan
            if(timesRan == nwarmup)   epsilon <<- exp(logEpsilonBar)
        },
        adapt_M = function() {
            ## logic follows Stan, but we use 1-based indexing
            in_adaptation_window <-
                (adaptWindow_counter > adapt_initBuffer) &
                (adaptWindow_counter <= nwarmup - adapt_termBuffer) &
                (adaptWindow_counter != nwarmup + 1)   ## seems redundant, but following Stan closely
            ##
            if(in_adaptation_window)   warmupSamples[adaptWindow_iter, 1:d] <<- state_sample$q
            end_adaptation_window <- (adaptWindow_counter == adapt_next_window) & (adaptWindow_counter != nwarmup + 1)  # ditto comment
            if(end_adaptation_window) {
                origM <- M
                for(i in 1:d) {
                    v <- var(warmupSamples[1:adaptWindow_iter, i])
                    vReg <- (adaptWindow_iter/(adaptWindow_iter+5))*v + 0.001*(5/(adaptWindow_iter+5))
                    M[i] <<- 1/vReg
                    sqrtM[i] <<- sqrt(M[i])
                }
                if(adapt_next_window == nwarmup - adapt_termBuffer) {
                    setSize(warmupSamples, 0, 0)   ## done, no further adaptation
                } else {
                    adaptWindow_size <<- adaptWindow_size * 2
                    adapt_next_window <<- adaptWindow_counter + adaptWindow_size
                    if(adapt_next_window != nwarmup - adapt_termBuffer) {
                        next_window_boundary <- adapt_next_window + 2*adaptWindow_size
                        if(next_window_boundary > nwarmup - adapt_termBuffer) {
                            adapt_next_window <<- nwarmup - adapt_termBuffer
                            adaptWindow_size <<- adapt_next_window - adaptWindow_counter
                        }
                    }
                    setSize(warmupSamples, adaptWindow_size, d, copy = FALSE, fillZeros = FALSE)
                }
                adaptWindow_iter <<- 1
                adaptWindow_counter <<- adaptWindow_counter + 1
                return(TRUE)
            }
            if(in_adaptation_window)   adaptWindow_iter <<- adaptWindow_iter + 1
            adaptWindow_counter <<- adaptWindow_counter + 1
            returnType(logical())
            return(FALSE)
        },
        before_chain = function(MCMCniter = double(), MCMCnburnin = double(), MCMCchain = double()) {
            if(MCMCchain == 1)   nwarmup <<- my_setWarmup$run(MCMCniter, MCMCnburnin, adaptive)
            if(adaptive) {
                if(nwarmup > 0) {
                    ## https://mc-stan.org/docs/2_23/reference-manual/hmc-algorithm-parameters.html#adaptation.figure
                    ## https://discourse.mc-stan.org/t/new-adaptive-warmup-proposal-looking-for-feedback/12039
                    ## https://colcarroll.github.io/hmc_tuning_talk/
                    ## approach follows Stan code
                    if(initBuffer + adaptWindow + termBuffer > nwarmup) {
                        if(messages & adaptive) print('  [Warning] Number of warmup iterations for NUTS_classic sampler ',
                                                      'is too small for even one cycle of standard adaptation. Using 15% ',
                                                      'for initial stepsize adaptation, 75% for mass matrix and stepsize ',
                                                      'adaptation, and 10% for final stepsize adaptation.')
                        adapt_initBuffer <<- round(nwarmup * 0.15)
                        adapt_termBuffer <<- round(nwarmup * 0.10)
                        adaptWindow_size <<- nwarmup - adapt_initBuffer - adapt_termBuffer
                    } else {
                        adaptWindow_size <<- adaptWindow
                        adapt_initBuffer <<- initBuffer
                        adapt_termBuffer <<- termBuffer
                        ## if there won't be room for the next window of doubled size, make the first and only window longer
                        if((nwarmup - (adapt_initBuffer + adaptWindow_size + adapt_termBuffer)) < 2*adaptWindow_size)
                            adaptWindow_size <<- nwarmup - (adapt_initBuffer + adapt_termBuffer)
                    }
                    # if(nwarmup < 20 & adaptive) if(messages) print("  [Warning] Number of warmup iteration for NUTS sampler is so small (",nwarmup,") that it might be useless.")
                    adapt_next_window <<- adapt_initBuffer + adaptWindow_size
                    adaptWindow_counter <<- 1
                    adaptWindow_iter <<- 1
                    Hbar <<- 0
                    logEpsilonBar <<- 0
                    stepsizeCounter <<- 0
                    setSize(warmupSamples, adaptWindow_size, d, fillZeros = FALSE)
                }
            }
        },
        after_chain = function() {
            if(messages) {
                if(numDivergences == 1)        print('  [Note] NUTS sampler (nodes: ', targetNodesToPrint, ') encountered ', numDivergences, ' divergent path.')
                if(numDivergences  > 1)        print('  [Note] NUTS sampler (nodes: ', targetNodesToPrint, ') encountered ', numDivergences, ' divergent paths.')
                if(numTimesMaxTreeDepth == 1)  print('  [Note] NUTS sampler (nodes: ', targetNodesToPrint, ') reached the maximum search tree depth ', numTimesMaxTreeDepth, ' time.')
                if(numTimesMaxTreeDepth  > 1)  print('  [Note] NUTS sampler (nodes: ', targetNodesToPrint, ') reached the maximum search tree depth ', numTimesMaxTreeDepth, ' times.')
                numDivergences <<- 0           ## reset counters for numDivergences and numTimesMaxTreeDepth,
                numTimesMaxTreeDepth <<- 0     ## even when using reset=FALSE to continue the same chain
            }
            if(warningInd > 0) {
                for(i in 1:warningInd) {
                    if(warningCodes[i,1] == 1) print('  [Warning] NUTS sampler (nodes: ', targetNodesToPrint, ') encountered a NaN value on MCMC iteration ', warningCodes[i,2], '.')
                    if(warningCodes[i,1] == 2) print('  [Warning] NUTS sampler (nodes: ', targetNodesToPrint, ') encountered acceptance prob = NaN in initEpsilon routine.')
                    if(warningCodes[i,1] == 3) print('  [Warning] NUTS sampler (nodes: ', targetNodesToPrint, ') encountered epsilon = NaN on MCMC iteration ', warningCodes[i,2], '.')
                }
                warningInd <<- 0               ## reset warningInd even when using reset=FALSE to continue the same chain
            }
        },
        reset = function() {
            timesRan       <<- 0
            epsilon        <<- epsilonOrig
            mu             <<- 0
            logEpsilonBar  <<- 0
            Hbar           <<- 0
            numDivergences <<- 0
            numTimesMaxTreeDepth <<- 0
            warningInd     <<- 0
            M              <<- Morig
            sqrtM          <<- sqrt(M)
            ## the adapt_* variables are initialized in before_chain()
            adaptWindow_size    <<- 0
            adapt_initBuffer    <<- 0
            adapt_termBuffer    <<- 0
            adapt_next_window   <<- 0
            adaptWindow_counter <<- 0
            adaptWindow_iter    <<- 0
            stepsizeCounter     <<- 0
        }
    ),
    buildDerivs = list(
        inverseTransformStoreCalculate = list(),
        calcLogProb = list(),
        gradient_aux = list()
    )
)


