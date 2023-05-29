


## #' Langevin Sampler
## #'
## #' The langevin sampler implements a special case of Hamiltonian Monte Carlo (HMC) sampling where only a single leapfrog step is taken on each sampling iteration, and the leapfrog step-size is adapted to match the scale of the posterior distribution (independently for each dimension being sampled).  The single leapfrog step is done by introducing auxiliary momentum variables, and using first-order derivatives to simulate Hamiltonian dynamics on this augmented paramter space (Neal, 2011).  Langevin sampling can operate on one or more continuous-valued posterior dimensions.  This sampling technique is also known as Langevin Monte Carlo (LMC), and the Metropolis-Adjusted Langevin Algorithm (MALA).
## #'
## #' @param model An uncompiled nimble model object on which the MCMC will operate.
## #' @param mvSaved A nimble \code{modelValues} object to be used to store MCMC samples.
## #' @param target A character vector of node names on which the sampler will operate.
## #' @param control A named list that controls the precise behavior of the sampler.  The default values for control list elements are specified in the setup code of the sampler.  A description of the possible control list elements appear in the details section.
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
## #' nimbleOptions(enableDerivs = TRUE)
## #' 
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



#' Hamiltonian Monte Carlo (HMC) Sampler
#'
#' The HMC sampler implements the No-U-Turn algorithm (NUTS; Hoffman and Gelman, 2014) for performing joint updates of multiple continuous-valued posterior dimensions.  This is done by introducing auxiliary momentum variables, and using first-order derivatives to simulate Hamiltonian dynamics on this augmented paramter space.  Internally, any posterior dimensions with bounded support are transformed, so sampling takes place on an unconstrained space.  In contrast to standard HMC (Neal, 2011), the NUTS algorithm removes the tuning parameters of the leapfrog step size and the number of leapfrog steps, thus providing a sampling algorithm that can be used without hand tuning or trial runs.
#'
#' @param model An uncompiled nimble model object on which the MCMC will operate.
#' @param mvSaved A nimble \code{modelValues} object to be used to store MCMC samples.
#' @param target A character vector of node names on which the sampler will operate.
#' @param control A named list that controls the precise behavior of the sampler.  The default values for control list elements are specified in the setup code of the sampler.  A description of the possible control list elements appear in the details section.
#'
#' @details
#'
#' The HMC sampler accepts the following control list elements:
#' 
#' \itemize{
#' \item messages.  A logical argument, specifying whether to print informative messages (default = TRUE)
#' \item numWarnings.  A numeric argument, specifying how many warnings messages to emit (for example, when NaN values are encountered).  See additional details below.  (default = 0)
#' \item gamma.  A positive numeric argument, specifying the degree of shrinkage used during the initial period of step-size adaptation. (default = 0.05)
#' \item initialEpsilon.  A positive numeric argument, specifying the initial step-size value. If not provided, an appropriate initial value is selected.
#' \item t0.  A non-negative numeric argument, where larger values stabilize (attenuate) the initial period of step-size adaptation. (default = 10)
#' \item kappa.  A numeric argument between zero and one, where smaller values give a higher weighting to more recent iterations during the initial period of step-size adaptation. (default = 0.75)
#' \item delta.  A numeric argument, specifying the target acceptance probability used during the initial period of step-size adaptation. (default = 0.65)
#' \item deltaMax.  A positive numeric argument, specifying the maximum allowable divergence from the Hamiltonian value. Paths which exceed this value are considered divergent, and will not proceed further. (default = 1000)
#' \item M.  A vector of positive real numbers, with length equal to the number of dimensions being sampled by HMC sampler.  Elements of M specify the diagonal elements of the diagonal mass matrix (or the metric) used for the auxiliary momentum variables in HMC sampling.  Sampling may be improved if the elements of M approximate the marginal inverse-variance (precision) the posterior dimensions.  (default: a vector of ones).
#' \item nwarmup.  The number of sampling iterations to adapt the leapfrog step-size.  This defaults to half the number of MCMC iterations, up to a maximum of 1000.
#' \item maxTreeDepth.  The maximum allowable depth of the binary leapfrog search tree for generating candidate transitions. (default = 10)
#' }
#'
#' NaN vales may be encountered in the course of the HMC leapfrog procedure.  In particular, when the stepsize (epsilon) is too large, the leapfrog procedure can step too far and arrive at an invalid region of parameter space, thus generating a NaN value in the likelihood evaluation or in the gradient calculation.  These situation are handled by the sampler by rejecting the NaN value, and reducing the stepsize.
#' 
#' @import nimble
#' 
#' @export
#'
#' @return A object of class `sampler_HMC`.
#' 
#' @aliases HMC hmc
#' 
#' @author Daniel Turek
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
#' conf$addSampler(target = c('b0', 'b1', 'sigma'), type = 'HMC')
#' 
#' Rmcmc <- buildMCMC(conf)
sampler_HMC <- nimbleFunction(
    name = 'sampler_HMC',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        printTimesRan  <- extractControlElement(control, 'printTimesRan',  FALSE)
        printEpsilon   <- extractControlElement(control, 'printEpsilon',   FALSE)
        printJ         <- extractControlElement(control, 'printJ',         FALSE)
        printM         <- extractControlElement(control, 'printM',         FALSE)
        messages       <- extractControlElement(control, 'messages',       getNimbleOption('verbose'))
        numWarnings    <- extractControlElement(control, 'numWarnings',    0)
        initialEpsilon <- extractControlElement(control, 'initialEpsilon', 0)
        gamma          <- extractControlElement(control, 'gamma',          0.05)
        t0             <- extractControlElement(control, 't0',             10)
        kappa          <- extractControlElement(control, 'kappa',          0.75)
        delta          <- extractControlElement(control, 'delta',          0.65)
        deltaMax       <- extractControlElement(control, 'deltaMax',       1000)
        M              <- extractControlElement(control, 'M',              -1)
        nwarmup        <- extractControlElement(control, 'nwarmup',        -1)
        maxTreeDepth   <- extractControlElement(control, 'maxTreeDepth',   10)
        ## node list generation
        targetNodes <- model$expandNodeNames(target)
        if(length(targetNodes) <= 0) stop('HMC sampler must operate on at least one node', call. = FALSE)
        targetNodesAsScalars <- model$expandNodeNames(targetNodes, returnScalarComponents = TRUE)
        targetNodesToPrint <- paste(targetNodes, collapse = ', ')
        if(nchar(targetNodesToPrint) > 100)   targetNodesToPrint <- paste0(substr(targetNodesToPrint, 1, 97), '...')
        calcNodes <- model$getDependencies(targetNodes)
        ## check for discrete nodes (early, before parameterTransform is specialized)
        if(any(model$isDiscrete(targetNodesAsScalars))) stop(paste0('HMC sampler cannot operate on discrete-valued nodes: ', paste0(targetNodesAsScalars[model$isDiscrete(targetNodesAsScalars)], collapse = ', ')))
        ## processing of bounds and transformations
        my_parameterTransform <- parameterTransform(model, targetNodesAsScalars)
        d <- my_parameterTransform$getTransformedLength()
        d2 <- max(d, 2) ## for pre-allocating vectors
        nimDerivs_wrt <- 1:d
        derivsInfo_return <- makeModelDerivsInfo(model, targetNodes, calcNodes)
        nimDerivs_updateNodes   <- derivsInfo_return$updateNodes
        nimDerivs_constantNodes <- derivsInfo_return$constantNodes
        ## numeric value generation
        timesRan <- 0;   epsilon <- 0;   mu <- 0;   logEpsilonBar <- 0;   Hbar <- 0
        q <- numeric(d2);   qL <- numeric(d2);   qR <- numeric(d2);   qDiff <- numeric(d2);   qNew <- numeric(d2)
        p <- numeric(d2);   pL <- numeric(d2);   pR <- numeric(d2);   p2 <- numeric(d2);      p3 <- numeric(d2)
        grad <- numeric(d2);   gradFirst <- numeric(d2);   gradSaveL <- numeric(d2);   gradSaveR <- numeric(d2)
        log2 <- log(2)
        warningCodes <- array(0, c(max(numWarnings,1), 2))
        warningInd <- 0
        nwarmupOrig <- nwarmup
        warmupIntervalLengths <- rep(0,7)            ## length 7 vector
        warmupIntervalsAdaptM <- rep(0,7)            ## length 7 vector
        warmupIntervalNumber <- 0
        warmupIntervalCount <- 0
        warmupSamples <- array(0, c(2,d2))           ## 2xd array
        if(length(M) == 1) { if(M == -1) M <- rep(1, d2) else M <- c(M, 1) }
        Morig <- M
        sqrtM <- sqrt(M)
        numDivergences <- 0
        numTimesMaxTreeDepth <- 0
        ## nested function and function list definitions
        qpNLDef <- nimbleList(q  = double(1), p  = double(1))
        btNLDef <- nimbleList(q1 = double(1), p1 = double(1), q2 = double(1), p2 = double(1), q3 = double(1), n = double(), s = double(), a = double(), na = double())
        ## checks
        if(!isTRUE(nimbleOptions('enableDerivs')))   stop('must enable NIMBLE derivatives, set nimbleOptions(enableDerivs = TRUE)', call. = FALSE)
        if(!isTRUE(model$modelDef[['buildDerivs']])) stop('must set buildDerivs = TRUE when building model',  call. = FALSE)
        if(initialEpsilon < 0) stop('HMC sampler initialEpsilon must be positive', call. = FALSE)
        if(!all(M > 0)) stop('HMC sampler M must contain all positive elements', call. = FALSE)
        if(d == 1) if(length(M) != 2) stop('length of HMC sampler M must match length of HMC target nodes', call. = FALSE)
        if(d  > 1) if(length(M) != d) stop('length of HMC sampler M must match length of HMC target nodes', call. = FALSE)
        if(maxTreeDepth < 1) stop('HMC maxTreeDepth must be at least one ', call. = FALSE)
    },
    run = function() {
        ## No-U-Turn Sampler with Dual Averaging, Algorithm 6 from Hoffman and Gelman (2014)
        if(timesRan == 0) {
            if(nwarmup == -1) stop('HMC nwarmup was not set correctly')
            ## reduce all pre-allocated vectors to correct size (d)
            q <<- q[1:d];   qL <<- qL[1:d];   qR <<- qR[1:d];   qDiff <<- qDiff[1:d];   qNew <<- qNew[1:d]
            p <<- p[1:d];   pL <<- pL[1:d];   pR <<- pR[1:d];   p2 <<- p2[1:d];           p3 <<- p3[1:d]
            grad <<- grad[1:d];   gradFirst <<- gradFirst[1:d];   gradSaveL <<- gradSaveL[1:d];   gradSaveR <<- gradSaveR[1:d]
            M <<- M[1:d];         sqrtM <<- sqrtM[1:d]
            if(initialEpsilon == 0) { initializeEpsilon()                 ## no initialEpsilon value was provided
                                  } else { epsilon <<- initialEpsilon }   ## user provided initialEpsilon
            mu <<- log(10*epsilon)
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
            if(v == -1) { btNL <- buildtree(qL, pL, logu, v, j, epsilon, qpLogH, 1)        ## first call: first = 1
                          qL <<- btNL$q1;   pL <<- btNL$p1
                      } else { btNL <- buildtree(qR, pR, logu, v, j, epsilon, qpLogH, 1)   ## first call: first = 1
                               qR <<- btNL$q2;   pR <<- btNL$p2 }
            if(btNL$s == 1)   if(runif(1) < btNL$n / n)   qNew <<- btNL$q3
            n <- n + btNL$n
            qDiff <<- qR - qL
            ##s <- btNL$s * nimStep(inprod(qDiff, pL)) * nimStep(inprod(qDiff, pR))                      ## this line replaced with the next,
            if(btNL$s == 0) s <- 0 else s <- nimStep(inprod(qDiff, pL)) * nimStep(inprod(qDiff, pR))     ## which acccounts for NaN's in btNL elements
            if(j >= maxTreeDepth) s <- 0
            if(printJ) {   if(j == 1) cat('j = ', j) else cat(', ', j)
                           cat('(');   if(v==1) cat('R') else cat('L');   cat(')')
                           if(s != 1) print(' ')   }
            if(j >= maxTreeDepth) { numTimesMaxTreeDepth <<- numTimesMaxTreeDepth + 1 }
            j <- j + 1
            checkInterrupt()
        }
        inverseTransformStoreCalculate(qNew)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        if(timesRan <= nwarmup)   adaptiveProcedure(btNL$a, btNL$na)
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
        gradient_aux = function(qArg = double(1)) {
            derivsOutput <- nimDerivs(inverseTransformStoreCalculate(qArg), order = 1, wrt = nimDerivs_wrt, model = model, updateNodes = nimDerivs_updateNodes, constantNodes = nimDerivs_constantNodes)
            returnType(double(1))
            return(derivsOutput$jacobian[1, 1:d])
        },
        gradient = function(qArg = double(1)) {
            derivsOutput <- nimDerivs(gradient_aux(qArg), order = 0, wrt = nimDerivs_wrt, model = model, updateNodes = nimDerivs_updateNodes, constantNodes = nimDerivs_constantNodes)
            grad <<- derivsOutput$value
        },
        leapfrog = function(qArg = double(1), pArg = double(1), eps = double(), first = double(), v = double()) {
            ## Algorithm 1 from Hoffman and Gelman (2014)
            if(first == 1) { gradient(qArg)     ## member data 'grad' is set in gradient() method
                         } else { if(v ==  1) grad <<- gradSaveR
                                  if(v == -1) grad <<- gradSaveL
                                  if(v ==  2) grad <<- gradSaveL }
            p2 <<- pArg + eps/2 * grad
            q2 <-  qArg + eps   * p2/M          ## see comments in drawMomentumValues method
            gradFirst <<- grad
            gradient(q2)                        ## member data 'grad' is set in gradient() method
            p3 <<- p2   + eps/2 * grad
            if(first == 1) { if(v ==  1) { gradSaveL <<- gradFirst;   gradSaveR <<- grad }
                             if(v == -1) { gradSaveR <<- gradFirst;   gradSaveL <<- grad }
                             if(v ==  2) { gradSaveL <<- gradFirst                       }
                         } else { if(v ==  1) gradSaveR <<- grad
                                  if(v == -1) gradSaveL <<- grad }
            if(warningInd < numWarnings) if(is.nan.vec(c(q2, p3))) { warningInd <<- warningInd + 1; warningCodes[warningInd,1] <<- 1; warningCodes[warningInd,2] <<- timesRan } ## message code 1: print('  [Warning] HMC sampler (nodes: ', targetNodesToPrint, ') encountered a NaN value in leapfrog routine, with timesRan = ', timesRan)
            returnType(qpNLDef());   return(qpNLDef$new(q = q2, p = p3))
        },
        initializeEpsilon = function() {
            ## Algorithm 4 from Hoffman and Gelman (2014)
            savedCalcNodeValues <- values(model, calcNodes)
            q <<- my_parameterTransform$transform(values(model, targetNodes))
            p <<- numeric(d)        ## keep, sets 'p' to size d on first iteration
            drawMomentumValues()    ## draws values for p
            epsilon <<- 1
            qpNL <- leapfrog(q, p, epsilon, 1, 2)            ## v = 2 is a special case for initializeEpsilon routine
            while(is.nan.vec(qpNL$q) | is.nan.vec(qpNL$p)) {              ## my addition
                ##if(numWarnings > 0) { print('  [Warning] HMC sampler (nodes: ', targetNodesToPrint, ') encountered NaN while initializing step-size; recommend better initial values')
                ##                      print('            reducing initial step-size'); numWarnings <<- numWarnings - 1 }
                epsilon <<- epsilon / 1000                                ## my addition
                qpNL <- leapfrog(q, p, epsilon, 0, 2)                     ## my addition
            }                                                             ## my addition
            qpLogH <- logH(q, p)
            a <- 2*nimStep(exp(logH(qpNL$q, qpNL$p) - qpLogH) - 0.5) - 1
            if(warningInd < numWarnings) if(is.nan(a)) { warningInd <<- warningInd + 1; warningCodes[warningInd,1] <<- 2; warningCodes[warningInd,2] <<- timesRan } ## message code 2: print('  [Warning] HMC sampler (nodes: ', targetNodesToPrint, ') caught acceptance prob = NaN in initializeEpsilon routine')
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
            Hbar <<- (1 - 1/(timesRan+t0)) * Hbar + 1/(timesRan+t0) * (delta - a/na)
            logEpsilon <- mu - sqrt(timesRan)/gamma * Hbar
            epsilon <<- exp(logEpsilon)
            timesRanToNegativeKappa <- timesRan^(-kappa)
            logEpsilonBar <<- timesRanToNegativeKappa * logEpsilon + (1 - timesRanToNegativeKappa) * logEpsilonBar
            if(timesRan == nwarmup)   epsilon <<- exp(logEpsilonBar)
            if(warningInd < numWarnings) if(is.nan(epsilon)) { warningInd <<- warningInd + 1; warningCodes[warningInd,1] <<- 3; warningCodes[warningInd,2] <<- timesRan } ## message code 3: print('  [Warning] HMC sampler (nodes: ', targetNodesToPrint, ') value of epsilon is NaN, with timesRan = ', timesRan)
            ## adapt M:
            if(warmupIntervalNumber < 8) {
                warmupIntervalCount <<- warmupIntervalCount + 1
                if(warmupIntervalCount > warmupIntervalLengths[warmupIntervalNumber]) stop('something went wrong in HMC warmup book-keeping')
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
                    }
                    warmupIntervalCount <<- 0
                    warmupIntervalNumber <<- warmupIntervalNumber + 1
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
                ##           if(numWarnings > 0) { print('  [Warning] HMC sampler (nodes: ', targetNodesToPrint, ') encountered a divergent path on iteration ', timesRan, ', with divergence = ', logu - qpLogH)
                ##                                 numWarnings <<- numWarnings - 1 } }
                a <- min(1, exp(qpLogH - logH0))
                if(is.nan.vec(q) | is.nan.vec(p)) { n <- 0; s <- 0; a <- 0 }     ## my addition
                return(btNLDef$new(q1 = q, p1 = p, q2 = q, p2 = p, q3 = q, n = n, s = s, a = a, na = 1))
            } else {        ## recursively build left and right subtrees
                btNL1 <- buildtree(qArg, pArg, logu, v, j-1, eps, logH0, 0)
                if(btNL1$s == 1) {
                    if(v == -1) { btNL2 <- buildtree(btNL1$q1, btNL1$p1, logu, v, j-1, eps, logH0, 0)   ## recursive calls: first = 0
                                  btNL1$q1 <- btNL2$q1;   btNL1$p1 <- btNL2$p1
                              } else {
                                  btNL2 <- buildtree(btNL1$q2, btNL1$p2, logu, v, j-1, eps, logH0, 0)   ## recursive calls: first = 0
                                  btNL1$q2 <- btNL2$q2;   btNL1$p2 <- btNL2$p2 }
                    nSum <- btNL1$n + btNL2$n
                    if(nSum > 0)   if(runif(1) < btNL2$n / nSum)   btNL1$q3 <- btNL2$q3
                    qDiff <<- btNL1$q2-btNL1$q1
                    btNL1$a  <- btNL1$a  + btNL2$a
                    btNL1$na <- btNL1$na + btNL2$na
                    btNL1$s  <- btNL2$s * nimStep(inprod(qDiff, btNL1$p1)) * nimStep(inprod(qDiff, btNL1$p2))
                    btNL1$n  <- nSum
                }
                return(btNL1)
            }
        },
        before_chain = function(MCMCniter = double(), MCMCnburnin = double(), MCMCchain = double()) {
            if(nwarmup == -1)   nwarmup <<- min( floor(MCMCniter/2), 1000 )
            if(MCMCchain == 1) {
                if(messages) print('  [Note] HMC sampler (nodes: ', targetNodesToPrint, ') is using ', nwarmup, ' warmup iterations.')
                if(nwarmup <  80) { print('  [Error] HMC sampler (nodes: ', targetNodesToPrint, ') requires a minimum of 80 warmup iterations.'); stop() }
                if(nwarmup < 200) print('  [Warning] A minimum of 200 warmup iterations is recommended for HMC sampler (nodes: ', targetNodesToPrint, ').')
                if(nwarmup > MCMCniter) print('  [Warning] Running fewer MCMC iterations than number of HMC warmup iterations (nodes: ', targetNodesToPrint, ').')
            }
            ## https://mc-stan.org/docs/2_23/reference-manual/hmc-algorithm-parameters.html#adaptation.figure
            ## https://discourse.mc-stan.org/t/new-adaptive-warmup-proposal-looking-for-feedback/12039
            ## https://colcarroll.github.io/hmc_tuning_talk/
            warmupBaseInterval <- floor(nwarmup/40)                             ## stan: 25
            if(warmupBaseInterval < 1)   stop('HMC warmupBaseInterval not set correctly')
            warmupIntervalLengths <<- warmupBaseInterval * c(3, 1, 2, 4, 8, 20, 2)    ## stan: 75 | 25 | 50 | 100 | 200 | 500 | 50
            warmupIntervalsAdaptM <<- c(0, 1, 1, 1, 1, 1, 0)
            setSize(warmupSamples, 20*warmupBaseInterval, d)
        },
        after_chain = function() {
            if(messages) {
                if(numDivergences == 1) print('  [Note] HMC sampler (nodes: ', targetNodesToPrint, ') encountered ', numDivergences, ' divergent path.')
                if(numDivergences  > 1) print('  [Note] HMC sampler (nodes: ', targetNodesToPrint, ') encountered ', numDivergences, ' divergent paths.')
                if(numTimesMaxTreeDepth == 1) print('  [Note] HMC sampler (nodes: ', targetNodesToPrint, ') reached the maximum search tree depth ', numTimesMaxTreeDepth, ' time.')
                if(numTimesMaxTreeDepth  > 1) print('  [Note] HMC sampler (nodes: ', targetNodesToPrint, ') reached the maximum search tree depth ', numTimesMaxTreeDepth, ' times.')
                numDivergences <<- 0           ## reset counters for numDivergences and numTimesMaxTreeDepth,
                numTimesMaxTreeDepth <<- 0     ## even when using reset=FALSE to continue the same chain
            }
            if(warningInd > 0) {
                for(i in 1:warningInd) {
                    if(warningCodes[i,1] == 1) print('  [Warning] HMC sampler (nodes: ', targetNodesToPrint, ') encountered a NaN value on MCMC iteration ', warningCodes[i,2], '.')
                    if(warningCodes[i,1] == 2) print('  [Warning] HMC sampler (nodes: ', targetNodesToPrint, ') encountered acceptance prob = NaN in initializeEpsilon routine.')
                    if(warningCodes[i,1] == 3) print('  [Warning] HMC sampler (nodes: ', targetNodesToPrint, ') encountered epsilon = NaN on MCMC iteration ', warningCodes[i,2], '.')
                }
                warningInd <<- 0               ## reset warningInd even when using reset=FALSE to continue the same chain
            }
        },
        reset = function() {
            timesRan       <<- 0
            epsilon        <<- 0
            mu             <<- 0
            logEpsilonBar  <<- 0
            Hbar           <<- 0
            numDivergences <<- 0
            numTimesMaxTreeDepth <<- 0
            warningInd     <<- 0
            nwarmup        <<- nwarmupOrig
            M              <<- Morig
            sqrtM          <<- sqrt(M)
            warmupIntervalNumber <<- 1
            warmupIntervalCount  <<- 0
        }
    ),
    buildDerivs = list(
        inverseTransformStoreCalculate = list(),
        gradient_aux = list()
    )
)


