#!/usr/bin/env Rscript

requirements <- c(
    'igraph',
    'coda',
    'R6',
    'nimble',
    'testthat')

##library(devtools)
##devtools::install_github('nimble-dev/nimble', ref = 'devel', subdir = 'packages/nimble')

for(package in requirements) {
    if(!suppressPackageStartupMessages(require(package, character.only = TRUE)))
        install.packages(package, repos = 'http://cran.us.r-project.org')
}

