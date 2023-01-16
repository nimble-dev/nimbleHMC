

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
    samples <- runMCMC(Cmcmc, 100000)
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
    expect_no_error(Rmcmc <- buildMCMC(conf))
    ##
    code <- nimbleCode({ x <- 3; y ~ T(dnorm(0, 1), 1, x) })
    Rmodel <- nimbleModel(code, inits = list(y = 2), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'HMC')
    expect_no_error(Rmcmc <- buildMCMC(conf))
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
    expect_no_error(Rmcmc <- buildMCMC(conf))
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
    expect_no_error(Rmcmc <- buildMCMC(conf))
})


## copied from 'Dirichlet-multinomial conjugacy' test in test-mcmc.R
test_that('HMC for Dirichlet-multinomial', {
    nimbleOptions(enableDerivs = TRUE)
    nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
    ##
    set.seed(0)
    n <- 100
    alpha <- c(10, 30, 15, 60, 1)
    K <- length(alpha)
    p <- c(.12, .24, .09, .54, .01)
    y <- rmulti(1, n, p)
    ##
    code <- nimbleCode({
        y[1:K] ~ dmulti(p[1:K], n);
        p[1:K] ~ ddirch(alpha[1:K]);
        for(i in 1:K) {
            alpha[i] ~ dgamma(.001, .001);
        }
    })
    constants <- list(K = K, n = n)
    inits <- list(p = rep(1/K, K), alpha = rep(K, K))
    data <- list(y = y)
    ##
    Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
    ##
    conf <- configureMCMC(Rmodel, nodes = NULL, monitors = 'p')
    addHMC(conf, 'alpha')
    addHMC(conf, 'p')
    Rmcmc <- buildMCMC(conf)
    ##
    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
    Cmcmc <- compiledList$mcmc
    ##
    set.seed(0)
    samples <- runMCMC(Cmcmc, niter = 20000, nburnin = 10000)
    ##
    expect_equal(as.numeric(apply(samples, 2, mean)), p, tol = .05)
})


## copied from 'block sampler on MVN node' test in test-mcmc.R
test_that('HMC on MVN node', {
    nimbleOptions(enableDerivs = TRUE)
    nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
    ##
    code <- nimbleCode({
        mu[1] <- 10
        mu[2] <- 20
        mu[3] <- 30
        x[1:3] ~ dmnorm(mu[1:3], prec = Q[1:3,1:3])
    })
    Q <- matrix(c(1.0,0.2,-1.0,0.2,4.04,1.6,-1.0,1.6,10.81), nrow=3)
    constants <- list(Q = Q)
    data <- list()
    inits <- list(x = c(10, 20, 30))
    ##
    Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
    Rmcmc <- buildHMC(Rmodel)
    ##
    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
    Cmcmc <- compiledList$mcmc
    ##
    set.seed(0)
    samples <- runMCMC(Cmcmc, niter = 20000, nburnin = 10000)
    ##
    expect_equal(as.numeric(apply(samples, 2, mean)), c(10,20,30), tol = .001)
    expect_equal(as.numeric(apply(samples, 2, var)), diag(solve(Q)), tol = .02)
})


## copied from 'test of conjugate Wishart' test in test-mcmc.R
test_that('HMC on conjugate Wishart', {
    nimbleOptions(enableDerivs = TRUE)
    nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
    ##
    set.seed(0)
    trueCor <- matrix(c(1, .3, .7, .3, 1, -0.2, .7, -0.2, 1), 3)
    covs <- c(3, 2, .5)
    trueCov <- diag(sqrt(covs)) %*% trueCor %*% diag(sqrt(covs))
    Omega <- solve(trueCov)
    n <- 20
    R <- diag(rep(1,3))
    mu <- 1:3
    Y <- mu + t(chol(trueCov)) %*% matrix(rnorm(3*n), ncol = n)
    M <- 3
    code <- nimbleCode({
        for(i in 1:n) {
            Y[i, 1:M] ~ dmnorm(mu[1:M], Omega[1:M,1:M])
        }
        Omega[1:M,1:M] ~ dwish(R[1:M,1:M], 4)
    })
    constants <- list(n = n, M = M, mu = mu, R = R)
    data <- list(Y = t(Y))
    inits <- list(Omega = diag(M))
    ##
    Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
    Rmcmc <- buildHMC(Rmodel)
    ##
    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
    Cmcmc <- compiledList$mcmc
    ##
    set.seed(0)
    samples <- runMCMC(Cmcmc, niter = 1000)
    ##
    newDf <- 4 + n
    newR <- R + tcrossprod(Y- mu)
    OmegaTrueMean <- newDf * solve(newR)
    wishRV <- array(0, c(M, M, 10000))
    for(i in 1:10000) {
        z <- t(chol(solve(newR))) %*% matrix(rnorm(3*newDf), ncol = newDf)
        wishRV[ , , i] <- tcrossprod(z)
    }
    OmegaSimTrueSDs <- apply(wishRV, c(1,2), sd)
    ##
    expect_equal(as.numeric(apply(samples, 2, mean)), as.numeric(OmegaTrueMean), tol = 0.2)
    expect_equal(as.numeric(apply(samples, 2, sd)), as.numeric(OmegaSimTrueSDs), tol = 0.04)
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




