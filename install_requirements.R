#!/usr/bin/env Rscript

requirements <- c(
    'igraph',
    'coda',
    'R6',
##    'nimble',    ## XXXXXX add this line back in
    'testthat')

for(package in requirements) {
    if(!suppressPackageStartupMessages(require(package, character.only = TRUE)))
        install.packages(package, repos = 'http://cran.us.r-project.org')
}

## XXXXXXXXXXXXXXXXXXX remove below
install.packages('devtools', repos = 'http://cran.us.r-project.org')
library(devtools)
devtools::install_github('nimble-dev/nimble', ref = 'ADoak_without_HMC', subdir = 'packages/nimble')
## XXXXXXXXXXXXXXXXXX remove until here

