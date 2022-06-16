#!/usr/bin/env Rscript

requirements <- c(
    'igraph',
    'coda',
    'R6',
##    'nimble',    ## XXXXXX add this line back in
    'testthat')

for(package in requirements) {
    install.packages(package)
}

## XXXXXXXXXXXXXXXXX remove below
install.packages('devtools')
library(devtools)
devtools::install_github('nimble-dev/nimble', ref = 'ADoak', subdir = 'packages/nimble')
## XXXXXXXXXXXXXXXX remove until here


