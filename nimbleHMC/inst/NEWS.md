#                            CHANGES IN VERSION 0.2.3 (December 2024)

- Minor bug fix in the adapatation routines of the `NUTS` and `NUTS_classic` samplers, which could leave model in an inconsistent state.

#                            CHANGES IN VERSION 0.2.1 (February 2024)

- Added checks within `NUTS` and `NUTS_classic` samplers to check validity of target nodes.

- Testing of correct results of HMC samplers on CAR model structures.

- `NUTS` and `NUTS_classic` samplers updated consistent with version 1.1.0 of `nimble`.

#                            CHANGES IN VERSION 0.2.0 (September 2023)

- Renamed `HMC` sampler to `NUTS_classic`, and added `NUTS` HMC sampler which follows the implementation of version 2.32.2 of Stan.

- Replaced `nwarmup` control argument to HMC samplers with the new `warmupMode` control argument.  This provides four distinct methods of choosing the number of HMC warmup iterations.

#                            CHANGES IN VERSION 0.1.1 (July 2023)

- Fix bug in `configureHMC` for sampler assignment to posterior predictive nodes.

#                            VERSION 0.1.0 (May 2023)

Initial release of `nimbleHMC` package.  Initial version provides:

- Hamiltonian Monte Carlo (HMC) sampler for use with `nimble` package MCMC engine.

- Convenience functions `addHMC`, `configureHMC`, `buildHMC` and `nimbleHMC`, for managing MCMC configuration objects with HMC samplers.
