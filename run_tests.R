#!/usr/bin/env Rscript

##testFiles <-
##    grep('test-.+\\.R$',
##         list.files('~/github/nimble/nimbleHMC/nimbleHMC/tests/testthat',
##                    full.names = TRUE),
##         value = TRUE)
## 
##for(f in testFiles) {
##    source(f)
##}

testFiles <-
    grep('test-.+\\.R$',
         list.files('nimbleHMC/tests/testthat',
                    full.names = TRUE),
         value = TRUE)


library(methods)
library(testthat)
library(nimbleHMC)

for(test in testFiles) {
    source(test)
}

##runner <- 'Rscript'
## 
### Run under exec_wait if sys package is installed, to support CTRL+C interrupts.
##if (require(sys)) {
##    custom_system2 <- sys::exec_wait
##    custom_shQuote <- function(x) x  # exec_wait doesn't like shQuote.
##} else {
##    cat('Missing suggested package sys, falling back to system2\n')
##    custom_system2 <- system2
##    custom_shQuote <- shQuote
##}
## 
##runTest <- function(test) {
##    cat('--------------------------------------------------------------------------------\n')
##    cat('TESTING', test, '\n')
##    if (runViaTestthat) {
##        name <- gsub('test-(.*)\\.R', '\\1', test)
##        script <- paste0('library(methods);',
##                         'library(testthat);',
##                         'library(nimble);',
##                         'tryCatch(test_package("nimble", "^', name, '$",',
##                         '                      reporter = ', reporter, '),',
##                         '  error = function(e) quit(status = 1))')
##        command <- c(runner, '-e', custom_shQuote(script))
##    } else command <- c(runner, file.path('packages', 'nimble', 'tests', 'testthat', test))
##    Sys.setenv(MAKEFLAGS = '-j1')  # Work around broken job pipe when GNU make is run under mclapply.
##    if (logToFile) {
##        logDir <- '/tmp/log/nimble'
##        dir.create(logDir, recursive = TRUE, showWarnings = FALSE)
##        stderr.log <- file.path(logDir, paste0('test-', name, '.stderr'))
##        stdout.log <- file.path(logDir, paste0('test-', name, '.stdout'))
##        if (custom_system2(command[1], tail(command, -1), stderr.log, stdout.log)) {
##            cat('\x1b[31mFAILED\x1b[0m', test, 'See', stderr.log, stdout.log, '\n')
##            return(TRUE)
##        }
##    } else {
##        if (custom_system2(command[1], tail(command, -1))) {
##            stop(paste('\x1b[31mFAILED\x1b[0m', test))
##        }
##    }
##    cat('\x1b[32mPASSED\x1b[0m', test, '\n')
##    return(FALSE)
##}
## 
##for (test in allTests) {
##    runTest(test)
##}



