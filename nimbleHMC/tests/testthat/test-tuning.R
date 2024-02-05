# Some tests of the tuning schedules, specifically:
# - setup of the sequence of windows for tuning M (variance or mass matrix)
# - setup when nwarmup is small
# - allowing nwarmup=0 to run with no self-tuning

nimbleOptions(enableDerivs=TRUE) # TRUE by default, but just in case

temporarilyAssignInGlobalEnv <- function(value, replace = FALSE) {
    name <- deparse(substitute(value))
    assign(name, value, envir = .GlobalEnv)
    if(!replace) {
        rmCommand <- substitute(remove(name, envir = .GlobalEnv))
        do.call('on.exit', list(rmCommand, add = TRUE), envir = parent.frame())
    }
}

capture_windows_NUTS <- function(sNUTS, niter, nwarmup) {
  sNUTS$before_chain(niter, 0, 1)
  update <- rep(FALSE, niter)
  AWI <- rep(0, niter)
  for(timesRan in 1:niter) {
    sNUTS$state_sample$q <- rnorm(sNUTS$d)
    if(timesRan <= nwarmup) {
      AWI[timesRan] <- sNUTS$adaptWindow_iter
      update[timesRan] <- sNUTS$adapt_M()
    }
  }
  list(update = update, AWI = AWI)
}

test_that('variance (mass) adaptation windows are set correctly', {
  # This test uses uncompiled only.
  # Adaptation windows are set up by different implementations of the
  # same logic in NUTS vs NUTS_classic.
  # In NUTS, at each iteration it is dynamically determined whether
  # the iteration is in an adaptation window and whether it marks the
  # end of a window.  In NUTS_classic, the equivalent set of windows
  # are precomputed once at the beginning. That makes NUTS_classic
  # easier to test. To test NUTS, we have the capture_windows_NUTS above
  # which sets a dummy random draw of states (state_sample$q) that will be
  # used by adapt_M, and then we call adapt_M solely to get its
  # dynamic windowing behavior.
  #
  # These different implementations are entirely due to our nimbleHMC coding history
  # within nimbleHMC and unrelated to how NUTS or NUTS_classic sampling works.
  #
  # Both behaviors follow what is in Hoffman and Gelman, with specifics that are
  # unstated in the paper (and/or modified by experience) determined from the
  # Stan source code.
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

  # First check windows with the default nwarmup=1000
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 1000))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  windows <- capture_windows_NUTS(sNUTS, 2000, 1000)
  AWI <- windows$AWI
  update <- windows$update
  expect_equal(which(update), c(100, 150, 250, 450, 950))
  expect_equal(AWI[1:1000], c(rep(1, 75),
                              1:25, 1:50, 1:100, 1:200, 1:500,
                              rep(1, 50)))

  # Now check with a small nwarmup value, where the window should be
  # split 75% / 15% / 10% (see output message)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 50))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  windows <- capture_windows_NUTS(sNUTS, 2000, 50)
  AWI <- windows$AWI
  update <- windows$update
  expect_equal(which(update), c(45))
  expect_equal(AWI[1:50], c(rep(1, 8),
                            1:37,
                            rep(1, 5)))

  # Now check with a larger number where the window doubling
  # and buffer scheme should fit exactly in the nwarmup
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 1700))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  windows <- capture_windows_NUTS(sNUTS, 1800, 1700)
  AWI <- windows$AWI
  update <- windows$update
  expect_equal(which(update), c(100, 150, 250, 450, 850, 1650))
  expect_equal(AWI[1:1700], c(rep(1, 75),
                              1:25, 1:50, 1:100, 1:200, 1:400, 1:800,
                              rep(1, 50)))

  # Now check when nwarmup is below a perfect-fit number, so
  # there will be one less window and the last window will be very long.
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 1699))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  windows <- capture_windows_NUTS(sNUTS, 1800, 1699)
  AWI <- windows$AWI
  update <- windows$update
  expect_equal(which(update), c(100, 150, 250, 450, 1649))
  expect_equal(AWI[1:1699], c(rep(1, 75),
                              1:25, 1:50, 1:100, 1:200, 1:(400+800-1),
                              rep(1, 50)))

  # Now check when nwarmup is above a perfect-fit number,
  # so the last window will simply be a bit longer.
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 1701))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  windows <- capture_windows_NUTS(sNUTS, 1800, 1701)
  AWI <- windows$AWI
  update <- windows$update
  expect_equal(which(update), c(100, 150, 250, 450, 850, 1651))
  expect_equal(AWI[1:1701], c(rep(1, 75),
                              1:25, 1:50, 1:100, 1:200, 1:400, 1:801,
                              rep(1, 50)))

  ## Now check that giving 1 <= nwarmup <= 19 throws an error
  ## but nwarmup == 0 is allowed.
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 19))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  expect_error(sNUTS$before_chain(2000, 0, 1))
  sNUTS$before_chain(2000, 0, 2)

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 0))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  sNUTS$before_chain(2000, 0, 2)
  Rmcmc$run(niter = 1)

  ###################################################
  ## Now repeat all the above cases with NUTS_classic
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 1000))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  sNUTS$before_chain(2000, 0, 1)
  expect_equal(sNUTS$warmupIntervalLengths,
               c(75, 25, 50, 100, 200, 500, 50))
  expect_equal(sNUTS$warmupIntervalsAdaptM,
               c(0, rep(1, length(sNUTS$warmupIntervalLengths)-2), 0))

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 50))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  sNUTS$before_chain(2000, 0, 1)
  expect_equal(sNUTS$warmupIntervalLengths,
               c(8, 37, 5))
  expect_equal(sNUTS$warmupIntervalsAdaptM,
               c(0, rep(1, length(sNUTS$warmupIntervalLengths)-2), 0))

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 1700))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  sNUTS$before_chain(1800, 0, 1)
  expect_equal(sNUTS$warmupIntervalLengths,
               c(75, 25, 50, 100, 200, 400, 800, 50))
  expect_equal(sNUTS$warmupIntervalsAdaptM,
               c(0, rep(1, length(sNUTS$warmupIntervalLengths)-2), 0))

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 1699))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  sNUTS$before_chain(1800, 0, 1)
  expect_equal(sNUTS$warmupIntervalLengths,
               c(75, 25, 50, 100, 200, 1199, 50))
  expect_equal(sNUTS$warmupIntervalsAdaptM,
               c(0, rep(1, length(sNUTS$warmupIntervalLengths)-2), 0))

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 1701))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  sNUTS$before_chain(1800, 0, 1)
  expect_equal(sNUTS$warmupIntervalLengths,
               c(75, 25, 50, 100, 200, 400, 801, 50))
  expect_equal(sNUTS$warmupIntervalsAdaptM,
               c(0, rep(1, length(sNUTS$warmupIntervalLengths)-2), 0))

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 19))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  expect_error(sNUTS$before_chain(2000, 0, 1))
  sNUTS$before_chain(2000, 0, 2)

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 19))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  expect_error(sNUTS$before_chain(2000, 0, 1))
  sNUTS$before_chain(2000, 0, 2)

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 0))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  sNUTS$before_chain(2000, 0, 2)
  Rmcmc$run(niter = 1)
})

test_that('epsilon and M adaptation vs not cases are handled correctly', {
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

  # temporarilyAssignInGlobalEnv(Rmodel)
  #A
  set.seed(2)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 1000, adaptive = FALSE, initializeEpsilon = TRUE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  Rmcmc$run(niter = 1)
  expect_true(sNUTS$epsilon == 1)

  #B
  set.seed(1)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 1000, adaptive = FALSE, initializeEpsilon = FALSE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  Rmcmc$run(niter = 1)
  expect_true(sNUTS$epsilon == 1)

  #C
  set.seed(1)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 1000, adaptive = FALSE, epsilon = 0.5, initializeEpsilon = FALSE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  Rmcmc$run(niter = 1)
  expect_true(sNUTS$epsilon == 0.5)

  #D
  set.seed(1)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 1000, adaptive = FALSE, epsilon = 0.5, initializeEpsilon = TRUE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  Rmcmc$run(niter = 1)
  expect_true(sNUTS$epsilon == 0.5)

  #E
  set.seed(1)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 1000, adaptive = TRUE, epsilon = 0.5, initializeEpsilon = FALSE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  Rmcmc$run(niter = 1)
  expect_false(sNUTS$epsilon == 0.5)

  #F
  set.seed(1)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 20, adaptive = TRUE, epsilon = 0.5, adaptEpsilon = FALSE, initializeEpsilon = FALSE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  expect_true(all(sNUTS$M == 1))
  Rmcmc$run(niter = 20)
  expect_true(sNUTS$epsilon == 0.5)
  expect_false(all(sNUTS$M == 1))

  #G
  set.seed(1)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 1000, adaptive = TRUE, epsilon = 0.5, adaptM = FALSE, initializeEpsilon = FALSE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  expect_true(all(sNUTS$M == 1))
  Rmcmc$run(niter = 1)
  expect_false(sNUTS$epsilon == 0.5)
  expect_true(all(sNUTS$M == 1))

  #A
  set.seed(2)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 1000, adaptive = FALSE, initializeEpsilon = TRUE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  Rmcmc$run(niter = 1)
  expect_true(sNUTS$epsilon == 1)

  #B
  set.seed(1)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 1000, adaptive = FALSE, initializeEpsilon = FALSE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  Rmcmc$run(niter = 1)
  expect_true(sNUTS$epsilon == 1)

  #C
  set.seed(1)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 1000, adaptive = FALSE, epsilon = 0.5, initializeEpsilon = FALSE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  Rmcmc$run(niter = 1)
  expect_true(sNUTS$epsilon == 0.5)

  #D
  set.seed(1)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 1000, adaptive = FALSE, epsilon = 0.5, initializeEpsilon = TRUE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  Rmcmc$run(niter = 1)
  expect_true(sNUTS$epsilon == 0.5)

  #E
  set.seed(1)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 1000, adaptive = TRUE, epsilon = 0.5, initializeEpsilon = FALSE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  Rmcmc$run(niter = 1)
  expect_false(sNUTS$epsilon == 0.5)

  #F
  set.seed(1)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 20, adaptive = TRUE, epsilon = 0.5, adaptEpsilon = FALSE, initializeEpsilon = FALSE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  expect_true(all(sNUTS$M == 1))
  Rmcmc$run(niter = 20)
  expect_true(sNUTS$epsilon == 0.5)
  expect_false(all(sNUTS$M == 1))

  #G##G
  set.seed(1)
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'iterations', warmup = 1000, adaptive = TRUE, epsilon = 0.5, adaptM = FALSE, initializeEpsilon = FALSE))
  Rmcmc <- buildMCMC(conf)
  sNUTS <- Rmcmc$samplerFunctions[[1]]
  expect_true(all(sNUTS$M == 1))
  Rmcmc$run(niter = 1)
  expect_false(sNUTS$epsilon == 0.5)
  expect_true(all(sNUTS$M == 1))

  # For NUTS_classic, we also wrote versions with compiling
  # These may turn out to be redundant with above.
  nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
  Rmodel2 <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
  
  # temporarilyAssignInGlobalEnv(Rmodel2)
  conf <- configureMCMC(Rmodel2, nodes = NULL)
  conf$addSampler('a', "NUTS_classic", control = list(warmupMode = 'default', adaptive = FALSE, initializeEpsilon = TRUE))
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel2)
  Cmcmc <- compileNimble(Rmcmc, project=Rmodel2)

  #A
  set.seed(1)
  Cmcmc$run(niter = 1)
  expect_true(Cmcmc$samplerFunctions[[1]]$epsilon == 1)

  #B
  Cmcmc$samplerFunctions[[1]]$epsilon <- 1
  Cmcmc$samplerFunctions[[1]]$initializeEpsilon <- FALSE
  set.seed(1)
  Cmcmc$run(niter = 1)
  expect_true(Cmcmc$samplerFunctions[[1]]$epsilon == 1)

  #C
  Cmcmc$samplerFunctions[[1]]$epsilon <- 0.5
  Cmcmc$samplerFunctions[[1]]$epsilonOrig <- 0.5
  set.seed(1)
  Cmcmc$run(niter = 1)
  expect_true(Cmcmc$samplerFunctions[[1]]$epsilon == 0.5)

  #D
  Cmcmc$samplerFunctions[[1]]$adaptive <- FALSE
  Cmcmc$samplerFunctions[[1]]$initializeEpsilon <- TRUE
  Cmcmc$samplerFunctions[[1]]$epsilon <- 0.5
  Cmcmc$samplerFunctions[[1]]$epsilonOrig <- 0.5
  set.seed(1)
  Cmcmc$run(niter = 1)
  expect_true(Cmcmc$samplerFunctions[[1]]$epsilon == 0.5)

  #E
  Cmcmc$samplerFunctions[[1]]$adaptive <- TRUE
  Cmcmc$samplerFunctions[[1]]$initializeEpsilon <- FALSE
  Cmcmc$samplerFunctions[[1]]$epsilon <- 0.5
  Cmcmc$samplerFunctions[[1]]$epsilonOrig <- 0.5
  set.seed(1)
  Cmcmc$run(niter = 40)
  expect_false(Cmcmc$samplerFunctions[[1]]$epsilon == 0.5)

  ###F
  Cmcmc$samplerFunctions[[1]]$adaptive <- TRUE
  Cmcmc$samplerFunctions[[1]]$initializeEpsilon <- FALSE
  Cmcmc$samplerFunctions[[1]]$epsilon <- 0.5
  Cmcmc$samplerFunctions[[1]]$epsilonOrig <- 0.5
  Cmcmc$samplerFunctions[[1]]$adaptEpsilon <- FALSE
  Cmcmc$samplerFunctions[[1]]$M <- rep(1, 3)
  Cmcmc$samplerFunctions[[1]]$sqrtM <- rep(1, 3)
  set.seed(1)
  expect_true(all(Cmcmc$samplerFunctions[[1]]$M == 1))
  Cmcmc$run(niter = 40)
  expect_true(Cmcmc$samplerFunctions[[1]]$epsilon == 0.5)
  expect_false(all(Cmcmc$samplerFunctions[[1]]$M == 1))
  ##
  ###G
  Cmcmc$samplerFunctions[[1]]$adaptive <- TRUE
  Cmcmc$samplerFunctions[[1]]$initializeEpsilon <- FALSE
  Cmcmc$samplerFunctions[[1]]$epsilon <- 0.5
  Cmcmc$samplerFunctions[[1]]$epsilonOrig <- 0.5
  Cmcmc$samplerFunctions[[1]]$adaptEpsilon <- TRUE
  Cmcmc$samplerFunctions[[1]]$adaptM <- FALSE
  Cmcmc$samplerFunctions[[1]]$M <- rep(1, 3)
  Cmcmc$samplerFunctions[[1]]$sqrtM <- rep(1, 3)
  set.seed(1)
  expect_true(all(Cmcmc$samplerFunctions[[1]]$M == 1))
  Cmcmc$run(niter = 40)
  expect_false(Cmcmc$samplerFunctions[[1]]$epsilon == 0.5)
  expect_true(all(Cmcmc$samplerFunctions[[1]]$M == 1))
})

test_that('burnin/warmup are handled correctly', {
    
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
  Cmodel <- compileNimble(Rmodel)

  ## 'default'
  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS")
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  expect_output(out <- runMCMC(Cmcmc, nburnin = 200, niter = 1000),
         regex = "using 200 warmup iterations.*all samples returned will be post-warmup")

  expect_output(out <- runMCMC(Cmcmc, niter = 1000),
         regex = "using 500 warmup iterations.*so the first half of the samples returned.*are from the warmup period")

  ## 'burnin'
  Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
  Cmodel <- compileNimble(Rmodel)

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'burnin'))
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  expect_output(out <- runMCMC(Cmcmc, niter = 1000),
         regex = "using 0 warmup iterations.*No adaptation is being done")

  Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
  Cmodel <- compileNimble(Rmodel)

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'burnin'))
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  expect_output(out <- runMCMC(Cmcmc, nburnin = 200, niter = 1000),
         regex = "using 200 warmup iterations.*all samples returned will be post-warmup")

  ## 'fraction'
  Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
  Cmodel <- compileNimble(Rmodel)

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'fraction', warmup = 0.25))
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  expect_output(out <- runMCMC(Cmcmc, niter = 1000),
                regex = "using 250 warmup iterations.*some of the samples returned will be collected during the warmup period")
  
  Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
  Cmodel <- compileNimble(Rmodel)

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'fraction', warmup = 0.25))
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  expect_output(out <- runMCMC(Cmcmc, nburnin = 100, niter = 1000),
                regex = "using 250 warmup iterations.*some of the samples returned will be collected during the warmup period")

  Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
  Cmodel <- compileNimble(Rmodel)

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'fraction', warmup = 0.25))
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  expect_output(out <- runMCMC(Cmcmc, nburnin = 250, niter = 1000),
                regex = "using 250 warmup iterations.*all samples returned will be post-warmup")

  ## 'iterations'
  Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
  Cmodel <- compileNimble(Rmodel)

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 500))
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  expect_output(out <- runMCMC(Cmcmc, niter = 1000),
                regex = "using 500 warmup iterations.*some of the samples returned will be collected during the warmup period")

  Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
  Cmodel <- compileNimble(Rmodel)

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(warmupMode = 'iterations', warmup = 500))
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  expect_output(out <- runMCMC(Cmcmc, nburnin = 500, niter = 1000),
                regex = "using 500 warmup iterations.*all samples returned will be post-warmup")

  ## no adaptation
  Rmodel <- nimbleModel(code, constants, data, inits, buildDerivs = TRUE)
  Cmodel <- compileNimble(Rmodel)

  conf <- configureMCMC(Rmodel, nodes = NULL)
  conf$addSampler('a', "NUTS", control = list(adaptive = FALSE))

  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel, resetFunctions = TRUE)

  expect_output(out <- runMCMC(Cmcmc, niter = 1000),
         regex = "has adaptation turned off")

})


