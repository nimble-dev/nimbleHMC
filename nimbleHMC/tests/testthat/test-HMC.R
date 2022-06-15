

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
    samples <- runMCMC(Cmcmc, 10000)
    ##
    expect_true(all(round(as.numeric(samples[1000,]), 5) == c(-0.11556, 0.88505, 2.89503)))
    expect_true(all(round(as.numeric(apply(samples, 2, mean)), 7) == c(0.4136721, 1.8417534, 3.2569954)))
    expect_true(all(round(as.numeric(apply(samples, 2, sd)), 7) == c(0.9268822, 1.2033593, 1.3179108)))
})


test_that('HMC sampler exact samples for different maxTreeDepths', {
    nimbleOptions(enableDerivs = TRUE)
    nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
    ##
    code <- nimbleCode({
        sigma ~ dunif(0, 100)
        a ~ dnorm(0, 0.01)
        y1 ~ dnorm(a, sd = sigma)
        y2 ~ dnorm(a, sd = 2*sigma)
    })
    constants <- list()
    data <- list(y1 = 1, y2 = 10)
    inits <- list(sigma = 1, a = 0)
    ##
    Rmodel1 <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
    conf1 <- configureMCMC(Rmodel1)
    conf1$removeSamplers(c('sigma', 'a'))
    conf1$addSampler(c('sigma', 'a'), type = 'HMC', control = list(nwarmup = 1000))
    Rmcmc1 <- buildMCMC(conf1)
    ##
    Rmodel2 <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
    conf2 <- configureMCMC(Rmodel2)
    conf2$removeSamplers(c('sigma', 'a'))
    conf2$addSampler(c('sigma', 'a'), type = 'HMC',
                     control = list(nwarmup = 1000, maxTreeDepth = 4))
    Rmcmc2 <- buildMCMC(conf2)
    ##
    Rmodel3 <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
    conf3 <- configureMCMC(Rmodel3)
    conf3$removeSamplers('sigma')
    conf3$addSampler('sigma', type = 'HMC', control = list(nwarmup = 1000))
    Rmcmc3 <- buildMCMC(conf3)
    ##
    compiledList <- compileNimble(list(model1=Rmodel1, mcmc1=Rmcmc1, model2=Rmodel2, mcmc2=Rmcmc2, model3=Rmodel3, mcmc3=Rmcmc3))
    Cmodel1 <- compiledList$model1; Cmcmc1 <- compiledList$mcmc1
    Cmodel2 <- compiledList$model2; Cmcmc2 <- compiledList$mcmc2
    Cmodel3 <- compiledList$model3; Cmcmc3 <- compiledList$mcmc3
    ##
    set.seed(0);   samples1 <- runMCMC(Cmcmc1, 10000)
    expect_true(all(round(as.numeric(samples1[10000,]),6) == c(-4.153855, 33.602766)))
    set.seed(0);   samples2 <- runMCMC(Cmcmc2, 10000)
    expect_true(all(round(as.numeric(samples2[10000,]),6) == c(-6.563781, 93.015175)))
    set.seed(0);   samples3 <- runMCMC(Cmcmc3, 10000)
    expect_true(all(round(as.numeric(samples3[10000,]),6) == c(-6.187742, 8.602518)))
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


test_that('HMC sampler reports correct number of divergences and max tree depths', {
    nimbleOptions(enableDerivs = TRUE)
    nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
    ##
    code <- nimbleCode({
        mu ~ dnorm(0, sd = 10)
        sigma ~ dunif(0, 10)
        p ~ dunif(0, 1)
        for(i in 1:3) {
            y[i] ~ dnorm(mu, sd = sigma)
        }
        yb ~ dbinom(size = 10, prob = p)
    })
    constants <- list()
    data <- list(y = c(10, 11, 19), yb = 7)
    inits <- list(mu = 0, sigma = 1, p = 0.5)
    Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
    ##
    conf <- configureMCMC(Rmodel)
    conf$addSampler(target = c('mu','sigma','p'), type = 'HMC',
                    control = list(nwarmup = 1000, maxTreeDepth = 5))
    Rmcmc <- buildMCMC(conf)
    ##
    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
    Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc
    ##
    set.seed(0)
    samples <- runMCMC(Cmcmc, 10000)
    ##
    expect_equal(Cmcmc$samplerFunctions$contentsList[[4]]$numDivergences,         9)
    expect_equal(Cmcmc$samplerFunctions$contentsList[[4]]$numTimesMaxTreeDepth, 699)
})




