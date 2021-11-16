
library(devtools)
library(nimble)

setwd('~/github/nimble/nimbleHMC')

document('nimbleHMC')

system('R CMD BUILD nimbleHMC')

check('nimbleHMC')

suppressMessages(try(remove.packages('nimbleHMC'), silent = TRUE))
tarFiles <- grep('\\.tar\\.gz', list.files(), value = TRUE)
lastTarFile <- tarFiles[length(tarFiles)]
message('installing package version ', gsub('\\.tar\\.gz$', '', lastTarFile))
system(paste0('R CMD install ', lastTarFile))

q('no')

library(nimbleHMC)

?tempFunction
tempFunction(3)

?myTempSampler

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
inits <- list(b0=0, b1=10, sigma=100)
Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$setSamplers()
####conf$addSampler(target = c('b0', 'b1', 'sigma'), type = 'RW_block')  ## XXXXXX MAKE 'HMC'
conf$addSampler(target=c('b0','b1','sigma'), type='myTempSampler', scalarComponents=TRUE)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)
