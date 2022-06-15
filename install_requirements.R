#!/usr/bin/env Rscript

## XXX REMOVE:
######library_install_path <- .libPaths()[1L]

## XXXXXX remove this:
###print('RUNNING INSTALL_REQUIREMENTS.R')
###print(paste0('library_install_path: ', library_install_path))
###print('BBBB')
##
###install.packages('igraph')
###a <- library(igraph, lib.loc = library_install_path)
###print('installed igraph:')
###print(a)
#####
###install.packages('coda')
###a <- library(coda, lib.loc = library_install_path)
###print('installed coda:')
###print(a)



requirements <- c(
    'igraph',
    'coda',
    'R6',
##    'nimble',    ## XXXXXX add this line back in
    'testthat')

for(package in requirements) {
    ##if(!suppressPackageStartupMessages(require(package, character.only = TRUE)))
    ##install.packages(package, repos = 'http://cran.us.r-project.org')
    ##install.packages(package, type = 'source')
    install.packages(package)
}

## XXXXXXXXXXXXXXXXX remove below
install.packages('devtools')
library(devtools)
devtools::install_github('nimble-dev/nimble', ref = 'ADoak', subdir = 'packages/nimble')
## XXXXXXXXXXXXXXXX remove until here




## XXXXXX remove this:
print('CCCCC')
###a <- library(igraph, lib.loc = library_install_path)
###a <- library(coda, lib.loc = library_install_path)
###a <- library(R6, lib.loc = library_install_path)
###a <- library(testthat, lib.loc = library_install_path)
a <- library(igraph)
a <- library(coda)
a <- library(R6)
a <- library(testthat)
a <- library(nimble)
print('loaded packages:')
print(a)


