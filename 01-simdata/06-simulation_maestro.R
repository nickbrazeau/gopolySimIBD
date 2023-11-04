## .................................................................................
## Purpose: Maestro for running polySimIBD Simulations
##
## Author: Nick Brazeau
##
## Date: 27 October, 2023
##
## Notes:
## .................................................................................
library(tidyverse)
squarecoords <- readRDS("01-simdata/00-simulation_setup/inputs/squarecoords.rds")
migmatdf <- readRDS("01-simdata/00-simulation_setup/inputs/migmat_framework.RDS")
#++++++++++++++++++++++++++++++++++++++++++
###   Part 0: Magic Numbers: immutable throughout simulations     ####
#++++++++++++++++++++++++++++++++++++++++++
reps <- 1:100 # number of simulation realizations to perform
tlim <- 25 # assume IBD to 25 generations for recent coalescent
# Miles et al. 2016 (PMC5052046) & Taylor et al. 2019 (PMC6707449) gives us a recombination rate by 7.4e-7 M/bp
# Aimee gets this number by taking the inverse of Mile's estimate of the CO recombination rate of 13.5 kb/cM
rho <- 7.4e-7

# approximate average of Pf3d7 Chromosome Lengths
pflen <- 1.664e6
# assuming single chromosome for ease
# assuming 1e3 loci
pos <- list(sort(sample(pflen, 1e3)))

# set previously
nDemes <- nrow(squarecoords)
modname <- migmatdf$modname

#++++++++++++++++++++++++++++++++++++++++++
### Part 1: 1st Level (Ind) Varying Parameters        ####
#++++++++++++++++++++++++++++++++++++++++++
lambdaCOI <-  readRDS("01-simdata/00-simulation_setup/inputs/optim_lambda.RDS")
mscale <- c(0, 0.25, 0.5, 0.75, 1) # mix of superinfection vs coinfection \

#++++++++++++++++++++++++++++++++++++++++++
### Part 2: 2nd Level (Spatial) Varying Parameters        ####
#++++++++++++++++++++++++++++++++++++++++++
baseNe <- c(10, 25, 50, 75, 100)

#++++++++++++++++++++++++++++++++++++++++++
### Part 3: Expand out Varying Parameters        ####
#++++++++++++++++++++++++++++++++++++++++++
simdf <- tidyr::expand_grid(
  nDemes,
  modname,
  tlim,
  rho,
  pos,
  reps,
  lambdaCOI,
  mscale,
  baseNe
)

#++++++++++++++++++++++++++++++++++++++++++
### Part 4: JOIN 2nd Level (Spatial) Varying Parameters        ####
#++++++++++++++++++++++++++++++++++++++++++
NeVary <- readRDS("01-simdata/00-simulation_setup/inputs/NeVary_migmat_framework.RDS")
MigMatVary <- readRDS("01-simdata/00-simulation_setup/inputs/MigVary_migmat_framework.RDS")
# bring together
simdf <- dplyr::full_join(simdf, NeVary, by = "modname", relationship = "many-to-many")
simdf <- dplyr::full_join(simdf, MigMatVary, by = "modname", relationship = "many-to-many")
