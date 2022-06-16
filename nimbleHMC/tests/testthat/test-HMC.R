

test_that('HMC sampler seems to work', {
    nimbleOptions(enableDerivs = TRUE)
    nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
    code <- nimbleCode({
        a[1] ~ dnorm(0, 1)
        a[2] ~ dnorm(a[1]+1, 1)
        a[3] ~ dnorm(a[2]+1, 1)
        d ~ dnorm(a[3], sd=2)
    })
    constants <- list()
    data <- list(d = 5)
    inits <- list(a = rep(0, 3))
    Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
    Rmodel$calculate()
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('a', 'HMC', control = list(nwarmup = 1000))
    Rmcmc <- buildMCMC(conf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    set.seed(0)
    samples <- runMCMC(Cmcmc, niter = 100000)
    ##
    expect_true(all(abs(as.numeric(apply(samples, 2, mean)) - c(0.4288181, 1.8582433, 3.2853841)) < 0.01))
    expect_true(all(abs(as.numeric(apply(samples, 2, sd)) - c(0.9248042, 1.1964343, 1.3098622)) < 0.01))
})


test_that('HMC sampler error messages for transformations with non-constant bounds', {
    nimbleOptions(enableDerivs = TRUE)
    nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
    ##
    code <- nimbleCode({ x ~ dexp(1); y ~ dunif(1, x) })
    Rmodel <- nimbleModel(code, inits = list(x = 10), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'HMC')
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    code <- nimbleCode({ x ~ dexp(1); y ~ dunif(x, 10) })
    Rmodel <- nimbleModel(code, inits = list(x = 1), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'HMC')
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    code <- nimbleCode({ x ~ dexp(1); y ~ T(dnorm(0, 1), 0, x) })
    Rmodel <- nimbleModel(code, inits = list(x = 1), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'HMC')
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    code <- nimbleCode({ x <- 2+c; y ~ T(dnorm(0, 1), 0, x) })
    Rmodel <- nimbleModel(code, constants = list(c = 1), inits = list(y = 1), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'HMC')
    expect_error(Rmcmc <- buildMCMC(conf), NA)   ## means: expect_no_error
    ##
    code <- nimbleCode({ x <- 3; y ~ T(dnorm(0, 1), 1, x) })
    Rmodel <- nimbleModel(code, inits = list(y = 2), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'HMC')
    expect_error(Rmcmc <- buildMCMC(conf), NA)   ## means: expect_no_error
    ##
    code <- nimbleCode({ x ~ dexp(1); y ~ T(dnorm(0, 1), 1, x) })
    Rmodel <- nimbleModel(code, inits = list(x = 3, y = 2), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'HMC')
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    code <- nimbleCode({ x ~ dexp(1); y ~ T(dnorm(0, 1), x, 10) })
    Rmodel <- nimbleModel(code, inits = list(x = 1, y = 2), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'HMC')
    expect_error(Rmcmc <- buildMCMC(conf))
})


test_that('HMC sampler error messages for invalid M mass matrix arguments', {
    nimbleOptions(enableDerivs = TRUE)
    nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
    ##
    code <- nimbleCode({
        for(i in 1:5)    x[i] ~ dnorm(0, 1)
    })
    Rmodel <- nimbleModel(code, inits = list(x = rep(0, 5)), buildDerivs = TRUE)
    ##
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('x[1]', 'HMC', M = 4)
    expect_error(Rmcmc <- buildMCMC(conf), NA)   ## means: expect_no_error
    ##
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('x[1]', 'HMC', M = 0)
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('x[1:3]', 'HMC', M = 4)
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('x[1:3]', 'HMC', M = c(1,2))
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('x[1:3]', 'HMC', M = c(1,0,2))
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('x[1:3]', 'HMC', M = c(1,2,3))
    expect_error(Rmcmc <- buildMCMC(conf), NA)   ## means: expect_no_error
})
 
 
####test_that('HMC sampler reports correct number of divergences and max tree depths', {
####    nimbleOptions(enableDerivs = TRUE)
####    nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
####    ##
####    code <- nimbleCode({
####        mu ~ dnorm(0, sd = 10)
####        sigma ~ dunif(0, 10)
####        p ~ dunif(0, 1)
####        for(i in 1:3) {
####            y[i] ~ dnorm(mu, sd = sigma)
####        }
####        yb ~ dbinom(size = 10, prob = p)
####    })
####    constants <- list()
####    data <- list(y = c(10, 11, 19), yb = 7)
####    inits <- list(mu = 0, sigma = 1, p = 0.5)
####    Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
####    ##
####    conf <- configureMCMC(Rmodel)
####    conf$addSampler(target = c('mu','sigma','p'), type = 'HMC',
####                    control = list(nwarmup = 1000, maxTreeDepth = 5))
####    Rmcmc <- buildMCMC(conf)
####    ##
####    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
####    Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc
####    ##
####    set.seed(0)
####    samples <- runMCMC(Cmcmc, 10000)
####    ##
####    expect_equal(Cmcmc$samplerFunctions$contentsList[[4]]$numDivergences,         9)
####    expect_equal(Cmcmc$samplerFunctions$contentsList[[4]]$numTimesMaxTreeDepth, 699)
####})




