#!/usr/bin/env Rscript

requirements <- c(
    'igraph',
    'coda',
    'R6',
##    'nimble',    ## XXXXXX add this line back in
    'testthat')

for(package in requirements) {
    if(!suppressPackageStartupMessages(require(package, character.only = TRUE)))
        ##install.packages(package, repos = 'http://cran.us.r-project.org')
        install.packages(package, type = 'source')
}

## XXXXXXXXXXXXXXXXXXX remove below
install.packages('devtools', type = 'source')
library(devtools)
devtools::install_github('nimble-dev/nimble', ref = 'ADoak_without_HMC', subdir = 'packages/nimble')
## XXXXXXXXXXXXXXXXXX remove until here

