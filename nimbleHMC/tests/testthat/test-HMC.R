nimbleOptions(enableDerivs = TRUE)

temporarilyAssignInGlobalEnv <- function(value, replace = FALSE) {
    name <- deparse(substitute(value))
    assign(name, value, envir = .GlobalEnv)
    if(!replace) {
        rmCommand <- substitute(remove(name, envir = .GlobalEnv))
        do.call('on.exit', list(rmCommand, add = TRUE), envir = parent.frame())
    }
}

test_that('HMC sampler seems to work', {
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

test_that('HMC sampler on subset of nodes', {
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
    conf <- configureMCMC(Rmodel, nodes = 'a[2]')
    conf$addSampler(c('a[1]','a[3]'), 'HMC', control = list(nwarmup = 1000))
    Rmcmc <- buildMCMC(conf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    set.seed(1)
    samples <- runMCMC(Cmcmc, 100000)
    ##
    expect_true(all(abs(as.numeric(apply(samples, 2, mean)) - c(0.4288181, 1.8582433, 3.2853841)) < 0.01))
    expect_true(all(abs(as.numeric(apply(samples, 2, sd)) - c(0.9248042, 1.1964343, 1.3098622)) < 0.01))
})

test_that('HMC sampler error messages for transformations with non-constant bounds', {
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
    expect_equal(as.numeric(apply(samples, 2, var)), diag(solve(Q)), tol = .03)
})


## copied from 'test of conjugate Wishart' test in test-mcmc.R
test_that('HMC on conjugate Wishart', {
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


test_that('HMC on LKJ', {
    R <- matrix(c(
        1, 0.9, .3, -.5, .1,
        0.9, 1, .15, -.3, .1,
        .3, .15, 1, .3, .1,
        -.5, -.3, .3, 1, .1,
        .1,.1,.1,.1, 1)
      , 5, 5)

    U <- chol(R)

    sds <- c(5, 4, 3, 2, 1)

    set.seed(1)
    Sigma <- diag(sds)%*%R%*%diag(sds)

    n <- 100
    p <- 5
    y <- t(t(chol(Sigma))%*%matrix(rnorm(p*n),p,n))

    uppertri_mult_diag <- nimbleFunction(
        run = function(mat = double(2), vec = double(1)) {
            returnType(double(2))
            p <- length(vec)
            out <- matrix(nrow = p, ncol = p, init = FALSE)
            for(i in 1:p)
                out[ , i] <- mat[ , i] * vec[i]
            return(out)
        },
        buildDerivs = list(run = list(ignore = 'i'))
    )
    temporarilyAssignInGlobalEnv(uppertri_mult_diag)

    code <- nimbleCode({
        for(i in 1:n)
            y[i, 1:p] ~ dmnorm(mu[1:p], cholesky = U[1:p, 1:p], prec_param = 0)
        U[1:p,1:p] <- uppertri_mult_diag(Ustar[1:p, 1:p], sds[1:p])
        Ustar[1:p,1:p] ~ dlkj_corr_cholesky(1.3, p)
    })
    m <- nimbleModel(code, constants = list(n = n, p = p, mu = rep(0, p)),
                     data = list(y = y), inits = list(sds = sds, Ustar = U),
                     buildDerivs = TRUE)
    cm <- compileNimble(m)

    conf <- configureMCMC(m, nodes = NULL)
    conf$addSampler('Ustar', 'HMC')
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project = m)

    out <- runMCMC(cmcmc, 10000)
    outSigma <- matrix(0, nrow(out), p*p)
    for(i in 1:nrow(outSigma))
        outSigma[i,] <- t(matrix(out[i,], p, p)) %*% matrix(out[i,],p,p)
    
    stan_means <- c(1.00000000, 0.87580832, 0.41032781, -0.56213296, 0.09006483, 0.87580832,
                    1.00000000, 0.18682787, -0.33699708, 0.12656145, 0.41032781, 0.18682787,
                    1.00000000, 0.11984278, 0.10919301, -0.56213296, -0.33699708, 0.11984278,
                    1.00000000, 0.10392069, 0.09006483, 0.12656145, 0.10919301, 0.10392069,
                    1.00000000)
    stan_sds <- c(0.000000e+00, 1.789045e-02, 6.244945e-02, 5.393811e-02, 7.928870e-02,
                  1.789045e-02, 0.000000e+00, 8.376820e-02, 7.448420e-02, 8.411652e-02,
                  6.244945e-02, 8.376820e-02, 8.600611e-17, 8.132228e-02, 9.242809e-02,
                  5.393811e-02, 7.448420e-02, 8.132228e-02, 8.711701e-17, 8.605078e-02,
                  7.928870e-02, 8.411652e-02, 9.242809e-02, 8.605078e-02, 1.227811e-16)
    
    nim_means_block <- apply(outSigma[1001:nrow(out), ], 2, mean)
    nim_sds_block <- apply(outSigma[1001:nrow(out), ], 2, sd)
    
    cols <- matrix(1:(p*p), p, p)
    cols <- cols[upper.tri(cols)]
    
    expect_lt(max(abs(stan_means[cols] - nim_means_block[cols])),  0.005)
    expect_lt(max(abs(stan_sds[cols] - nim_sds_block[cols])), 0.005)
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


test_that('testing HMC configuration functions', {
    ##
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
    ##
    ## addHMC <- function(conf, nodes = character(), control = list(), replace = FALSE, print = TRUE) {}
    ##
    conf <- configureMCMC(Rmodel)
    addHMC(conf, nodes = 'a[1]')
    samp <- conf$getSamplers('a[1]')
    expect_equal(length(samp), 2)
    expect_equal(samp[[1]]$name, 'conjugate_dnorm_dnorm_additive')
    expect_equal(samp[[2]]$name, 'HMC')
    ##
    addHMC(conf, nodes = c('a[1]', 'a[3]'), replace = TRUE)
    expect_equal(length(conf$getSamplers()), 2)
    samp <- conf$getSamplers('a[1]')
    expect_equal(length(samp), 1)
    expect_equal(samp[[1]]$name, 'HMC')
    samp <- conf$getSamplers('a[3]')
    expect_equal(length(samp), 1)
    expect_equal(samp[[1]]$name, 'HMC')
    ##
    addHMC(conf, nodes = 'a', replace = TRUE)
    expect_equal(length(conf$getSamplers()), 1)
    expect_equal(length(conf$getSamplers('a[1]')), 1)
    expect_equal(length(conf$getSamplers('a[2]')), 1)
    expect_equal(length(conf$getSamplers('a[3]')), 1)
    ##
    ## configureHMC <- function(model, nodes = character(), control = list(), print = TRUE, ...) {}
    ##
    conf <- configureHMC(Rmodel)
    samp <- conf$getSamplers()
    expect_equal(length(samp), 1)
    expect_equal(samp[[1]]$name, 'HMC')
    expect_equal(samp[[1]]$target, c('a[1]', 'a[2]', 'a[3]'))
    ##
    conf <- configureHMC(Rmodel, nodes = c('a[1]', 'a[3]'))
    samp <- conf$getSamplers()
    expect_equal(length(samp), 1)
    expect_equal(samp[[1]]$name, 'HMC')
    expect_equal(samp[[1]]$target, c('a[1]', 'a[3]'))
    ##
    ## buildHMC <- function(model, nodes = character(), control = list(), print = TRUE, ...) {}
    ##
    Rmcmc <- buildHMC(Rmodel)
    samplers <- Rmcmc$samplerFunctions$contentsList
    expect_equal(length(samplers), 1)
    expect_equal(as.character(class(samplers[[1]])), 'sampler_HMC')
    expect_equal(samplers[[1]]$targetNodes, c('a[1]', 'a[2]', 'a[3]'))
    ##
    Rmcmc <- buildHMC(Rmodel, nodes = c('a[1]', 'a[3]'))
    samplers <- Rmcmc$samplerFunctions$contentsList
    expect_equal(length(samplers), 1)
    expect_equal(as.character(class(samplers[[1]])), 'sampler_HMC')
    expect_equal(samplers[[1]]$targetNodes, c('a[1]', 'a[3]'))
    ##
})


test_that('error trap discrete latent nodes', {
    code <- nimbleCode({
        y ~ dnorm(z, 1)
        z ~ dpois(mu)
        y2 ~ dnorm(z2, 1)
        z2 ~ dbinom(p, 5)
        mu ~ dnorm(0, 1)
        p ~ dunif(0, 1)
    })
    data <- list(y = 0.1, y2 = 0.3)
    inits <- list(z = 2, z2 = 3, mu = 0, p = .3)
    Rmodel <- nimbleModel(code, data = data, inits = inits, buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    addHMC(conf, nodes = 'mu')
    buildMCMC(conf)
    expect_error(addHMC(conf, nodes = 'z'), 'HMC sampler cannot be applied')
    expect_error(addHMC(conf, nodes = 'z2'), 'HMC sampler cannot be applied')
    expect_error(addHMC(conf, nodes = c('mu','z')), 'HMC sampler cannot be applied')
})


test_that('correctly assign samplers for discrete and continuous nodes', {
    code <- nimbleCode({
        mu ~ dnorm(0, 1)
        z ~ dpois(mu)
        y ~ dnorm(z, 1)
        p ~ dbeta(1, 1)
        z2 ~ dbinom(p, 5)
        y2 ~ dnorm(z2, 1)
    })
    data <- list(y = 0.1, y2 = 0.3)
    inits <- list(z = 2, z2 = 3, mu = 0, p = .3)
    Rmodel <- nimbleModel(code, data = data, inits = inits, buildDerivs = TRUE)
    ##
    conf <- configureHMC(Rmodel)
    expect_equal(length(conf$samplerConfs), 3)
    expect_identical(conf$samplerConfs[[1]]$target, c('mu', 'p'))
    expect_identical(conf$samplerConfs[[1]]$name, 'HMC')
    expect_identical(conf$samplerConfs[[2]]$target, 'z')
    expect_identical(conf$samplerConfs[[2]]$name, 'slice')
    expect_identical(conf$samplerConfs[[3]]$target, 'z2')
    expect_identical(conf$samplerConfs[[3]]$name, 'slice')
    ##
    conf <- configureHMC(Rmodel, nodes = c('z', 'mu'))
    expect_equal(length(conf$samplerConfs), 2)
    expect_identical(conf$samplerConfs[[1]]$target, 'mu')
    expect_identical(conf$samplerConfs[[1]]$name, 'HMC')
    expect_identical(conf$samplerConfs[[2]]$target, 'z')
    expect_identical(conf$samplerConfs[[2]]$name, 'slice')
    ##
    Rmcmc <- buildHMC(Rmodel)
    expect_equal(length(Rmcmc$samplerFunctions$contentsList), 3)
    expect_identical(Rmcmc$samplerFunctions$contentsList[[1]]$targetNodes, c('mu', 'p'))
    expect_identical(as.character(class(Rmcmc$samplerFunctions$contentsList[[1]])), 'sampler_HMC')
    expect_identical(Rmcmc$samplerFunctions$contentsList[[2]]$target, 'z')
    expect_identical(as.character(class(Rmcmc$samplerFunctions$contentsList[[2]])), 'sampler_slice')
    expect_identical(Rmcmc$samplerFunctions$contentsList[[3]]$target, 'z2')
    expect_identical(as.character(class(Rmcmc$samplerFunctions$contentsList[[3]])), 'sampler_slice')
    ##
    Rmcmc <- buildHMC(Rmodel, nodes = c('z', 'mu'))
    expect_equal(length(Rmcmc$samplerFunctions$contentsList), 2)
    expect_identical(Rmcmc$samplerFunctions$contentsList[[1]]$targetNodes, 'mu')
    expect_identical(as.character(class(Rmcmc$samplerFunctions$contentsList[[1]])), 'sampler_HMC')
    expect_identical(Rmcmc$samplerFunctions$contentsList[[2]]$target, 'z')
    expect_identical(as.character(class(Rmcmc$samplerFunctions$contentsList[[2]])), 'sampler_slice')
})


