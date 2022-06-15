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

## XXXXXXXXXXXXXXXXX remove below
install.packages('devtools', type = 'source')
library(devtools)
devtools::install_github('nimble-dev/nimble', ref = 'ADoak', subdir = 'packages/nimble')
## XXXXXXXXXXXXXXXX remove until here

