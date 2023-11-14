## .................................................................................
## Purpose: Maestro for running polySimIBD Simulations
##
## Author: Nick Brazeau
##
## Date: 27 October, 2023
##
## Notes:
## .................................................................................
set.seed(48)
library(tidyverse)
squarecoords <- readRDS("01-simdata/sim_params/squarecoords.rds")
migmatdf <- readRDS("01-simdata/sim_params/migmat_framework.RDS")
#++++++++++++++++++++++++++++++++++++++++++
###   Part 0: Magic Numbers: immutable throughout simulations     ####
#++++++++++++++++++++++++++++++++++++++++++
rep <- 1:100 # number of simulation realizations to perform
tlim <- 10 # assume IBD to 10 generations for recent coalescent
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
lambdaCOI <-  readRDS("01-simdata/sim_params/optim_lambda.RDS")
mscale <- c(0, 0.25, 0.5, 0.75, 1) # mix of superinfection vs coinfection

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
  rep,
  lambdaCOI,
  mscale,
  baseNe
)

#++++++++++++++++++++++++++++++++++++++++++
### Part 4: JOIN 2nd Level (Spatial) Varying Parameters        ####
#++++++++++++++++++++++++++++++++++++++++++
NeVary <- readRDS("01-simdata/sim_params/NeVary_migmat_framework.RDS")
MigMatVary <- readRDS("01-simdata/sim_params/MigVary_migmat_framework.RDS")
# bring together
simdf <- dplyr::full_join(simdf, NeVary, by = "modname", relationship = "many-to-many")
simdf <- dplyr::full_join(simdf, MigMatVary, by = "modname", relationship = "many-to-many")

#......................
# now read in migmat
#......................
getmigmat <- function(modname) {
  pth <- switch(modname,
                "torus" = "01-simdata/sim_params/migmat_surface_grids/torus.RDS",
                "dexter" = "01-simdata/sim_params/migmat_surface_grids/dexter.RDS",
                "bound" = "01-simdata/sim_params/migmat_surface_grids/bound.RDS"
         )
  return(readRDS(pth))
}

simdf <- simdf %>%
  dplyr::mutate(migmat = purrr::map(modname, getmigmat))



#++++++++++++++++++++++++++++++++++++++++++
### Part 5: Misc Pieces for Run       ####
#++++++++++++++++++++++++++++++++++++++++++
#......................
# get seeds
#......................
simdf$randseedkey <- floor(runif(n = nrow(simdf), min = 100, max = 1e6))
while ( length(unique(simdf$randseedkey)) != nrow(simdf)  ) {
  rwsfix <- which( duplicated(simdf$randseedkey) )
  simdf$randseedkey[rwsfix] <- floor(runif(n = length(rwsfix), min = 100, max = 1e6))
}

#......................
# down sample for IBD
#......................
simdf$dwnsmplnum <- 5


#++++++++++++++++++++++++++++++++++++++++++
### Part 6: Split for Inputs       ####
#++++++++++++++++++++++++++++++++++++++++++
# save out for posterity
saveRDS(simdf, "01-simdata/simulation_maestro.RDS")

#++++++++++++++++++++++++++++++++++++++++++
### Addendum 10/14/2023: Subset for DISCent cases    ####
#++++++++++++++++++++++++++++++++++++++++++
#TODO expand out in future
simdfOG <- simdf
simdf <- simdfOG %>%
  dplyr::filter(mscale == 0 & lambdaCOI < 1)

#......................
# now split
#......................
simdflist <- split(simdf, f = 1:nrow(simdf))
simdflistnm <- as.character( purrr::map_int(simdflist, "randseedkey") )
# now save out
dir.create("01-simdata/sim_inputs", recursive = T)
mapply(function(x,y) {
  saveRDS(x, file = paste0("01-simdata/sim_inputs/", "simswf_input_key-", y, ".RDS"))
  return(0)
}, x = simdflist, y = simdflistnm)



#++++++++++++++++++++++++++++++++++++++++++
### Part 7: Maestro TSV Map    ####
#++++++++++++++++++++++++++++++++++++++++++
maestro_snake_map <- tibble::tibble(
  input = paste0("01-simdata/sim_inputs/", "simswf_input_key", simdflistnm, ".RDS"),
  outpath = paste0("01-simdata/swfsim_results/", "simswf_result_key", simdflistnm, ".RDS")
)
readr::write_tsv(maestro_snake_map, "01-simdata/polySimIBD_simulationrun_maestro.tsv")


