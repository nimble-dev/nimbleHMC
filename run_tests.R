
testFiles <-
    grep('test-.+\\.R$',
         list.files('~/github/nimble/nimbleHMC/nimbleHMC/tests/testthat',
                    full.names = TRUE),
         value = TRUE)

for(f in testFiles) {
    source(f)
}





