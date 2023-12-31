## .................................................................................
## Purpose: Setup a NeVary Multiplier for Respective Surfaces
##
## Author: Nick Brazeau
##
## Date: 27 October, 2023
##
## Notes: Type 1 and Type 2 will have Ne Varying Sims
## .................................................................................
library(tidyverse)

#............................................................
# read in and manip
#...........................................................
migmatdf <- readRDS("01-simdata/sim_params/migmat_framework.RDS")
# replicate migmat type 1 and type 2
migmatdf <- dplyr::bind_rows(migmatdf, migmatdf[1:2,])
migmatdf$modnameNe <- NA
migmatdf$modnameNe[1:3] <- migmatdf$modname[1:3]
migmatdf$modnameNe[4] <- paste0(migmatdf$modname[4], "-NeVary")
migmatdf$modnameNe[5] <- paste0(migmatdf$modname[5], "-NeVary")
#......................
# set up Ne Vary column
#......................
migmatdf$NeVaryMult <- NA
# get number of inds
nInds <- nrow(migmatdf$migmat[[1]])
# vary by e^0.5*x, where x is 0:4
expmult <- exp(0.5 * 0:4)
expmult <- rep(expmult, sqrt(nInds))
migmatdf$NeVaryMult[1:3] <- list(rep(1, nInds))
migmatdf$NeVaryMult[4:5] <- list(expmult)

#............................................................
# out
#...........................................................
migmatdf <- migmatdf %>%
  dplyr::select(c("modnameNe", "modname", "NeVaryMult"))
saveRDS(migmatdf,
        "01-simdata/sim_params/NeVary_migmat_framework.RDS")
