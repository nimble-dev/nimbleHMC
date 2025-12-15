nimbleOptions(enableDerivs = TRUE) # TRUE by default, but just in case

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
    for(type in c('NUTS', 'NUTS_classic')) {
        cat(paste0('testing ', type, ' sampler\n'))
        Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
        Rmodel$calculate()
        conf <- configureMCMC(Rmodel, nodes = NULL)
        conf$addSampler('a', type, control = list(nwarmup = 1000))
        Rmcmc <- buildMCMC(conf)
        Cmodel <- compileNimble(Rmodel)
        Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
        set.seed(0)
        samples <- runMCMC(Cmcmc, 100000)
        ##
        expect_true(all(abs(as.numeric(apply(samples, 2, mean)) - c(0.4288181, 1.8582433, 3.2853841)) < 0.017))
        expect_true(all(abs(as.numeric(apply(samples, 2, sd)) - c(0.9248042, 1.1964343, 1.3098622)) < 0.011))
    }
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
    for(type in c('NUTS', 'NUTS_classic')) {
        cat(paste0('testing ', type, ' sampler\n'))
        Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
        Rmodel$calculate()
        conf <- configureMCMC(Rmodel, nodes = 'a[2]')
        conf$addSampler(c('a[1]','a[3]'), type = type, control = list(nwarmup = 1000))
        Rmcmc <- buildMCMC(conf)
        Cmodel <- compileNimble(Rmodel)
        Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
        set.seed(1)
        samples <- runMCMC(Cmcmc, 100000)
        expect_true(all(abs(as.numeric(apply(samples, 2, mean)) - c(0.4288181, 1.8582433, 3.2853841)) < 0.015))
        expect_true(all(abs(as.numeric(apply(samples, 2, sd)) - c(0.9248042, 1.1964343, 1.3098622)) < 0.013))
    }
})

test_that('HMC sampler error messages for transformations with non-constant bounds', {
    code <- nimbleCode({ x ~ dexp(1); y ~ dunif(1, x) })
    Rmodel <- nimbleModel(code, inits = list(x = 10), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'NUTS')
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    code <- nimbleCode({ x ~ dexp(1); y ~ dunif(x, 10) })
    Rmodel <- nimbleModel(code, inits = list(x = 1), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'NUTS')
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    code <- nimbleCode({ x ~ dexp(1); y ~ T(dnorm(0, 1), 0, x) })
    Rmodel <- nimbleModel(code, inits = list(x = 1), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'NUTS')
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    code <- nimbleCode({ x <- 2+c; y ~ T(dnorm(0, 1), 0, x) })
    Rmodel <- nimbleModel(code, constants = list(c = 1), inits = list(y = 1), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'NUTS')
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    code <- nimbleCode({ x <- 3; y ~ T(dnorm(0, 1), 1, x) })
    Rmodel <- nimbleModel(code, inits = list(y = 2), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'NUTS')
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    code <- nimbleCode({ x ~ dexp(1); y ~ T(dnorm(0, 1), 1, x) })
    Rmodel <- nimbleModel(code, inits = list(x = 3, y = 2), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'NUTS')
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    code <- nimbleCode({ x ~ dexp(1); y ~ T(dnorm(0, 1), x, 10) })
    Rmodel <- nimbleModel(code, inits = list(x = 1, y = 2), buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('y', 'NUTS')
    expect_error(Rmcmc <- buildMCMC(conf))
})

test_that('hmc_checkTarget catches all invalid cases', {
    code <- nimbleCode({
        x[1]   ~ dbern(0.5)
        x[2]   ~ dbin(size = 4, prob = 0.5)
        x[3]   ~ dcat(prob = p[1:3])
        x[4]   ~ dpois(2)
        x[5:7] ~ dmulti(prob = p[1:3], size = 3)
        x[8]   ~ T(dnorm(0, 1), 0,  )
        x[9]   ~ T(dnorm(0, 1),  , 2)
        x[10]  ~ T(dnorm(0, 1), 0, 2)
        ##
        a[1] ~ dnorm(0, 1)
        b[1] ~ T(dnorm(a[1], 1), 0, 2)
        ##
        a[2] ~ dnorm(0, 1)
        b[2] ~ dconstraint(a[2] > 0)
        ##
        a[3] ~ dnorm(0, 1)
        b[3] ~ dinterval(a[3], 0)
    })
    constants <- list(p = rep(1/3,3))
    inits <- list(x = rep(1, 10), a = rep(1,  3))
    data <- list(b = rep(1, 3))
    Rmodel <- nimbleModel(code, constants, data, inits)
    conf <- configureMCMC(Rmodel, nodes = NULL, print = FALSE)
    ##
    for(node in Rmodel$expandNodeNames(c('x', 'a'))) {
        conf$setSamplers()
        conf$addSampler(target = node, type = 'NUTS_classic')
        expect_error(buildMCMC(conf))
        conf$setSamplers()
        conf$addSampler(target = node, type = 'NUTS')
        expect_error(buildMCMC(conf))
    }
})

test_that('hmc_checkTarget catches non-AD support for custom distributions', {
    ddistNoAD <- nimbleFunction(
        run = function(x = double(0), log = integer(0, default = 0)) {
            returnType(double(0)); return(1)
        }
    )
    code <- nimbleCode({
        x ~ ddistNoAD()
        y ~ dnorm(x, 1)
    })
    Rmodel <- nimbleModel(code, data = list(y=0), inits = list(x=0), buildDerivs = TRUE)
    conf <- configureHMC(Rmodel)
    expect_error(buildMCMC(conf))
    ##
    ddistAD <- nimbleFunction(
        run = function(x = double(0), log = integer(0, default = 0)) {
            returnType(double(0)); return(1)
        },
        buildDerivs = TRUE
    )
    code <- nimbleCode({
        x ~ ddistAD()
        y ~ dnorm(x, 1)
    })
    Rmodel <- nimbleModel(code, data = list(y=0), inits = list(x=0), buildDerivs = TRUE)
    conf <- configureHMC(Rmodel)
    expect_no_error(buildMCMC(conf))
})

test_that('HMC sampler error messages for invalid M mass matrix arguments', {
    code <- nimbleCode({
        for(i in 1:5)    x[i] ~ dnorm(0, 1)
    })
    Rmodel <- nimbleModel(code, inits = list(x = rep(0, 5)), buildDerivs = TRUE)
    ##
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('x[1]', 'NUTS', M = 4)
    expect_no_error(Rmcmc <- buildMCMC(conf))
    ##
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('x[1]', 'NUTS', M = 0)
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('x[1:3]', 'NUTS', M = 4)
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('x[1:3]', 'NUTS', M = c(1,2))
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('x[1:3]', 'NUTS', M = c(1,0,2))
    expect_error(Rmcmc <- buildMCMC(conf))
    ##
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('x[1:3]', 'NUTS', M = c(1,2,3))
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
    for(type in c('NUTS', 'NUTS_classic')) {
        cat(paste0('testing ', type, ' sampler\n'))
        Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
        conf <- configureMCMC(Rmodel, nodes = NULL, monitors = 'p')
        addHMC(conf, type = type, 'alpha')
        addHMC(conf, type = type, 'p')
        Rmcmc <- buildMCMC(conf)
        compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
        Cmcmc <- compiledList$mcmc
        set.seed(0)
        samples <- runMCMC(Cmcmc, niter = 20000, nburnin = 10000)
        expect_equal(as.numeric(apply(samples, 2, mean)), y/n, tol = .02)
    }
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
    cat(paste0('testing NUTS_classic sampler\n'))
    Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    addHMC(conf, 'x', type = 'NUTS_classic')
    Rmcmc <- buildMCMC(conf)
    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
    Cmcmc <- compiledList$mcmc
    set.seed(0)
    samples <- runMCMC(Cmcmc, niter = 20000, nburnin = 10000)
    expect_equal(as.numeric(apply(samples, 2, mean)), c(10,20,30), tol = .001)
    expect_equal(as.numeric(apply(samples, 2, var)), diag(solve(Q)), tol = .03)
    ##
    cat(paste0('testing NUTS sampler\n'))
    Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, nodes = NULL)
    addHMC(conf, 'x', type = 'NUTS')
    Rmcmc <- buildMCMC(conf)
    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
    Cmcmc <- compiledList$mcmc
    set.seed(0)
    samples <- runMCMC(Cmcmc, niter = 20000, nburnin = 10000)
    expect_equal(as.numeric(apply(samples, 2, mean)), c(10,20,30), tol = .001)
    expect_equal(as.numeric(apply(samples, 2, var)), diag(solve(Q)), tol = .065)
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
    samples <- runMCMC(Cmcmc, niter = 5000, nburnin = 2000)
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
    expect_equal(as.numeric(apply(samples, 2, mean)), as.numeric(OmegaTrueMean), tol = 0.1)
    expect_equal(as.numeric(apply(samples, 2, sd)), as.numeric(OmegaSimTrueSDs), tol = 0.1)
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
    ##
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
    conf$addSampler('Ustar', 'NUTS')
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project = m)
    ##
    out <- runMCMC(cmcmc, 10000)
    outSigma <- matrix(0, nrow(out), p*p)
    for(i in 1:nrow(outSigma))
        outSigma[i,] <- t(matrix(out[i,], p, p)) %*% matrix(out[i,],p,p)
    ##
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
    ##
    nim_means_block <- apply(outSigma[1001:nrow(out), ], 2, mean)
    nim_sds_block <- apply(outSigma[1001:nrow(out), ], 2, sd)
    cols <- matrix(1:(p*p), p, p)
    cols <- cols[upper.tri(cols)]
    expect_lt(max(abs(stan_means[cols] - nim_means_block[cols])),  0.005)
    expect_lt(max(abs(stan_sds[cols] - nim_sds_block[cols])), 0.005)
})
 
test_that('testing HMC configuration functions', {
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
    conf <- configureMCMC(Rmodel)
    addHMC(conf, target = 'a[1]')
    samp <- conf$getSamplers('a[1]')
    expect_equal(length(samp), 2)
    expect_equal(samp[[1]]$name, 'conjugate_dnorm_dnorm_additive')
    expect_equal(samp[[2]]$name, 'NUTS')
    ##
    addHMC(conf, target = c('a[1]', 'a[3]'), replace = TRUE)
    expect_equal(length(conf$getSamplers()), 2)
    samp <- conf$getSamplers('a[1]')
    expect_equal(length(samp), 1)
    expect_equal(samp[[1]]$name, 'NUTS')
    samp <- conf$getSamplers('a[3]')
    expect_equal(length(samp), 1)
    expect_equal(samp[[1]]$name, 'NUTS')
    ##
    addHMC(conf, target = 'a', replace = TRUE)
    expect_equal(length(conf$getSamplers()), 1)
    expect_equal(length(conf$getSamplers('a[1]')), 1)
    expect_equal(length(conf$getSamplers('a[2]')), 1)
    expect_equal(length(conf$getSamplers('a[3]')), 1)
    ##
    conf <- configureHMC(Rmodel)
    samp <- conf$getSamplers()
    expect_equal(length(samp), 1)
    expect_equal(samp[[1]]$name, 'NUTS')
    expect_equal(samp[[1]]$target, c('a[1]', 'a[2]', 'a[3]'))
    ##
    conf <- configureHMC(Rmodel, nodes = c('a[1]', 'a[3]'))
    samp <- conf$getSamplers()
    expect_equal(length(samp), 1)
    expect_equal(samp[[1]]$name, 'NUTS')
    expect_equal(samp[[1]]$target, c('a[1]', 'a[3]'))
    ##
    Rmcmc <- buildHMC(Rmodel)
    samplers <- Rmcmc$samplerFunctions$contentsList
    expect_equal(length(samplers), 1)
    expect_equal(as.character(class(samplers[[1]])), 'sampler_NUTS')
    expect_equal(samplers[[1]]$targetNodes, c('a[1]', 'a[2]', 'a[3]'))
    ##
    Rmcmc <- buildHMC(Rmodel, nodes = c('a[1]', 'a[3]'))
    samplers <- Rmcmc$samplerFunctions$contentsList
    expect_equal(length(samplers), 1)
    expect_equal(as.character(class(samplers[[1]])), 'sampler_NUTS')
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
    addHMC(conf, target = 'mu')
    buildMCMC(conf)
    expect_error(addHMC(conf, target = 'z'), 'HMC sampler cannot be applied')
    expect_error(addHMC(conf, target = 'z2'), 'HMC sampler cannot be applied')
    expect_error(addHMC(conf, target = c('mu','z')), 'HMC sampler cannot be applied')
    expect_error(addHMC(conf, target = 'mu', type = 'wrong'))
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
    expect_identical(conf$samplerConfs[[1]]$name, 'NUTS')
    expect_identical(conf$samplerConfs[[2]]$target, 'z')
    expect_identical(conf$samplerConfs[[2]]$name, 'slice')
    expect_identical(conf$samplerConfs[[3]]$target, 'z2')
    expect_identical(conf$samplerConfs[[3]]$name, 'slice')
    ##
    conf <- configureHMC(Rmodel, nodes = c('z', 'mu'))
    expect_equal(length(conf$samplerConfs), 2)
    expect_identical(conf$samplerConfs[[1]]$target, 'mu')
    expect_identical(conf$samplerConfs[[1]]$name, 'NUTS')
    expect_identical(conf$samplerConfs[[2]]$target, 'z')
    expect_identical(conf$samplerConfs[[2]]$name, 'slice')
    ##
    Rmcmc <- buildHMC(Rmodel)
    expect_equal(length(Rmcmc$samplerFunctions$contentsList), 3)
    expect_identical(Rmcmc$samplerFunctions$contentsList[[1]]$targetNodes, c('mu', 'p'))
    expect_identical(as.character(class(Rmcmc$samplerFunctions$contentsList[[1]])), 'sampler_NUTS')
    expect_identical(Rmcmc$samplerFunctions$contentsList[[2]]$target, 'z')
    expect_identical(as.character(class(Rmcmc$samplerFunctions$contentsList[[2]])), 'sampler_slice')
    expect_identical(Rmcmc$samplerFunctions$contentsList[[3]]$target, 'z2')
    expect_identical(as.character(class(Rmcmc$samplerFunctions$contentsList[[3]])), 'sampler_slice')
    ##
    Rmcmc <- buildHMC(Rmodel, nodes = c('z', 'mu'))
    expect_equal(length(Rmcmc$samplerFunctions$contentsList), 2)
    expect_identical(Rmcmc$samplerFunctions$contentsList[[1]]$targetNodes, 'mu')
    expect_identical(as.character(class(Rmcmc$samplerFunctions$contentsList[[1]])), 'sampler_NUTS')
    expect_identical(Rmcmc$samplerFunctions$contentsList[[2]]$target, 'z')
    expect_identical(as.character(class(Rmcmc$samplerFunctions$contentsList[[2]])), 'sampler_slice')
})

test_that('configureHMC correctly assign samplers for posterior-predictive nodes', {
    code <- nimbleCode({
        mu ~ dnorm(0,1)
        sd ~ dunif(0, 10)
        y ~ dnorm(mu, sd = sd)
        pp ~ dnorm(mu, sd = sd)
    })
    constants <- list()
    data <- list(y = 0)
    inits <- list(mu = 0, sd = 1, pp = 0)
    Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
    conf <- configureHMC(Rmodel)
    ##
    expect_true(length(conf$samplerConfs) == 2)
    expect_true(conf$samplerConfs[[1]]$name == 'NUTS')
    expect_identical(conf$samplerConfs[[1]]$target, c('mu', 'sd'))
    expect_true(conf$samplerConfs[[2]]$name == 'posterior_predictive')
    expect_identical(conf$samplerConfs[[2]]$target, 'pp')
})

test_that('HMC runs with various non-differentiable constructs', {
    code <- nimbleCode({
        for(i in 1:3) {
            w[i] ~ dconstraint(theta[i] > 0)
            cens[i] ~ dinterval(t[i], c[i])
            z[i] ~ dcat(p[1:2])
            t[i] ~ dweib(r, 1)
            y[i] ~ dnorm(theta[i], 1)
            theta[i] ~ dnorm(theta0, 1)
        }
        p[1] ~ dunif(0,1)
        p[2] <- 1-p[1]
        r ~ dunif(0,5)
        theta0 ~ dnorm(0,1)
    })
    m <- nimbleModel(code, data = list(w = rep(1,3), t = c(NA,NA,0.5), cens = c(1,1,0), y = rnorm(3), z = c(1,1,2)),
                     constants = list(c=rep(1.5,3)), inits = list(t=c(2,3,NA), theta = runif(3,0,3), r = 1),
                     buildDerivs = TRUE)
    conf <- configureMCMC(m)
    conf$removeSampler('theta0')
    conf$removeSampler('r')
    conf$removeSampler('p[1]')
    conf$addSampler(c('theta0','r','p[1]'), 'NUTS')
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    out <- runMCMC(cmcmc, niter = 100)
    expect_identical(nrow(out), 100L)    
})

test_that('HMC results for CAR match non-HMC', {
    set.seed(1)
    code <- nimbleCode({
        S[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau)
        tau ~ dunif(0, 5)
        for(i in 1:N)
            mu[i] <- S[i]
        for(i in 1:N) {
            log(lambda[i]) <- mu[i]
            Y[i] ~ dpois(lambda[i])
        }
    })
    ##
    constants <- list(N = 6,
                      num = c(1,2,2,2,2,1),
                      adj = c(2, 1,3, 2,4, 3,5, 4,6, 5),
                      weights = rep(1, 10),
                      L = 10)
    data <- list(Y = c(1,0,2,1,4,3))
    inits <- list(tau = 1, S = c(0,0,0,0,0,0))
    ##
    Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
    conf <- configureMCMC(Rmodel, monitors = c('tau','S'))
    Rmcmc <- buildMCMC(conf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    out <- runMCMC(Cmcmc, niter = 505000, nburnin = 5000, thin = 500)
    ##
    Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
    conf <- configureHMC(Rmodel, monitors = c('tau','S'))
    Rmcmc <- buildMCMC(conf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    outHMC <- runMCMC(Cmcmc, niter = 22000, nburnin=2000, thin=20)
    ##
    expect_equal(apply(out[,1:6],2,mean), apply(outHMC[,1:6],2,mean), tolerance = .06)
    expect_equal(mean(out[,7]),mean(outHMC[,7]), tolerance = .15)
    expect_equal(apply(out,2,quantile,c(.1,.9)), apply(outHMC,2,quantile,c(.1,.9)), tolerance = 0.15)
})


test_that('HMC results for mixture model match non-HMC', {
    set.seed(1)
    code <- nimbleCode({
        for(i in 1:n) {
            y[i] ~ dnorm(mu[k[i]],1)
            k[i] ~ dcat(p[1:K])
        }
        for(i in 1:K)
            mu[i] ~ dnorm(mu0,1)
        p[1:K] ~ ddirch(alpha[1:K])
        mu0 ~ dflat()
    })
    ##
    n <- 500
    K <- 3
    constants <- list(n=n, K=K, alpha = rep(1,K))
    mu <- c(0,2,4)
    data <- list(y=sample(c(rnorm(50,mu[1],.35), rnorm(250,mu[2],.35), rnorm(200,mu[3],.35)), n, replace=FALSE))
    inits <- list(k=sample(1:K,n,replace=T),mu=c(-1,2,6),mu0=1)
    ##
    m <- nimbleModel(code, constants = constants, data = data,
                     inits = inits, buildDerivs = TRUE)
    conf <- configureMCMC(m, monitors=c('mu','mu0','p'))
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc)
    out <- runMCMC(cmcmc, niter=50000, nburnin=10000, thin=40)
    ##
    m <- nimbleModel(code, constants = constants, data = data,
                     inits = inits, buildDerivs = TRUE)
    conf <- configureMCMC(m, nodes=c('k'), monitors=c('mu','mu0','p'))
    conf$addSampler(c('mu0','mu','p'), 'NUTS')
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc)
    outHMC <- runMCMC(cmcmc, niter=22000, nburnin=2000, thin=20)
    ## deal with label switching
    sorter <- function(row) {
        ord <- order(row[1:3])
        return(c(row[1:3][ord], row[4], row[5:7][ord]))
    }
    out <- t(apply(out, 1, sorter))
    outHMC <- t(apply(outHMC, 1, sorter))
    ##
    expect_equal(apply(out,2,mean), apply(outHMC,2,mean), tolerance = 0.1)
    expect_equal(apply(out,2,quantile,c(.1,.9)), apply(outHMC,2,quantile,c(.1,.9)), tolerance = 0.1)
})


