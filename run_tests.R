#!/usr/bin/env Rscript

library(methods)
library(testthat)
library(nimble)
library(nimbleHMC)

testFiles <-
    grep('test-.+\\.R$',
         list.files('nimbleHMC/tests/testthat', full.names = TRUE),
         value = TRUE)

for(test in testFiles) {
    source(test)
}

