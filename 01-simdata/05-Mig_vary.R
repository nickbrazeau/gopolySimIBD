## .................................................................................
## Purpose: Setup a Migration Multiplier for Respective Surfaces
##
## Author: Nick Brazeau
##
## Date: 27 October, 2023
##
## Notes: Will want to vary how much time is spent at home
## .................................................................................
library(tidyverse)

#............................................................
# read in and manip
#...........................................................
migmatdf <- readRDS("01-simdata/sim_params/migmat_framework.RDS") %>%
  dplyr::select(-c("migmat"))
# I want to migration to vary such that home is preferred 99% ... 50%
MigVaryFct <- c(99, 95, 90, 75, 50)
# expand out
migmatdf <- tidyr::expand_grid(migmatdf, MigVaryFct)

#............................................................
# out
#...........................................................
migmatdf <- migmatdf %>%
  dplyr::select(c("modname", "MigVaryFct"))
saveRDS(migmatdf,
        "01-simdata/sim_params/MigVary_migmat_framework.RDS")
