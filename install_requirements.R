#!/usr/bin/env Rscript

## XXXXXX remove this:
print('RUNNING INSTALL_REQUIREMENTS.R')
install.packages('igraph')
a <- library(igraph)
print('installed igraph:')
print(a)



##### XXXXXXXXX uncomment below:
##requirements <- c(
##    'igraph',
##    'coda',
##    'R6',
####    'nimble',    ## XXXXXX add this line back in
##    'testthat')

####### XXXXXXXXX uncomment below:
##for(package in requirements) {
##    if(!suppressPackageStartupMessages(require(package, character.only = TRUE)))
##        ##install.packages(package, repos = 'http://cran.us.r-project.org')
##        install.packages(package, type = 'source')
##}

####### XXXXXXXXX uncomment below:
#### XXXXXXXXXXXXXXXXX remove below
##install.packages('devtools', type = 'source')
##library(devtools)
##devtools::install_github('nimble-dev/nimble', ref = 'ADoak', subdir = 'packages/nimble')
#### XXXXXXXXXXXXXXXX remove until here

