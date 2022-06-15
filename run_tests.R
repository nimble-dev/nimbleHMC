#!/usr/bin/env Rscript

### XXXX remove this:
print('RUNNING RUN_TESTS.R')

library(methods)
library(testthat)
library(nimble)
library(nimbleHMC)

testFiles <-
    grep('test-.+\\.R$',
         list.files('nimbleHMC/tests/testthat', full.names = TRUE),
         value = TRUE)

### XXXX uncomment below:
####for(test in testFiles) {
####    cat('===========================================================\n')
####    cat(paste0('Running test-', gsub('.*test-', '', test), '\n'))
####    cat('===========================================================\n')
####    source(test)
####}

