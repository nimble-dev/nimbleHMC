#!/usr/bin/env Rscript

requirements <- c(
    'igraph',
    'coda',
    'R6',
    'nimble',
    'testthat')

for(package in requirements) {
    install.packages(package)
}

