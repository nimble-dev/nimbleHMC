
test_that('tempSampler works', {
    ##
    code <- nimbleCode({
        b0 ~ dnorm(0, 0.001)
        b1 ~ dnorm(0, 0.001)
        sigma ~ dunif(0, 10000)
        for(i in 1:N) {
            mu[i] <- b0 + b1 * x[i]
            y[i] ~ dnorm(mu[i], sd = sigma)
        }
    })
    N <- 10
    constants <- list(N = N, x = 1:N)
    data <- list(y = 1:N)
    inits <- list(b0 = 1, b1 = 0.1, sigma = 1)
    Rmodel <- nimbleModel(code, constants, data, inits)
    lp <- Rmodel$calculate()
    expect_equal(round(lp,4), -138.5709)
    ##
    conf <- configureMCMC(Rmodel)
    conf$setSamplers()
    conf$addSampler(target = c('b0', 'b1', 'sigma'), type = 'RW_block')
    conf$addSampler(target = c('b0'), type = 'myTempSampler')
    conf$addSampler(target = c('b1'), type = 'myTempSampler')
    samp <- conf$getSamplers()
    expect_equal(length(samp), 3)
    expect_equal(samp[[1]]$name, 'RW_block')
    expect_equal(samp[[2]]$name, 'myTempSampler')
    expect_equal(samp[[3]]$name, 'myTempSampler')
    ##
    Rmcmc <- buildMCMC(conf)
    expect_true(class(Rmcmc) == 'MCMC')
    print('myTempSampler passed tests')
})

