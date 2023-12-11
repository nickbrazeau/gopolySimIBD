## .................................................................................
## Purpose: Wrapper script for running polySimIBD::sim_swf from command line
##
## Author: Nick Brazeau
## .................................................................................
# Temporarily suppress warning from RNG and "expect" val to be numeric for strict snakemake bash mode
defaultwarnings <- getOption("warn")
options(warn = -1)


#++++++++++++++++++++++++++++++++++++++++++
### dependencies     ####
#++++++++++++++++++++++++++++++++++++++++++
deps <- c("polySimIBD", "goodegg", "optparse", "dplyr", "tibble", "magrittr")
deps <- !sapply(deps, function(x){x %in% installed.packages()[,1]} ) # note this a named vector

# catch polySimIBD remote
if(deps["polySimIBD"]) {
  if (!"remotes" %in% installed.packages()[,1]){
    install.packages("remotes")
  }
  remotes::install_github("nickbrazeau/polySimIBD")
  deps <- deps[names(deps) != "polySimIBD"]
}

# catch goodegg remote
if(deps["goodegg"]) {
  if (!"remotes" %in% installed.packages()[,1]){
    install.packages("remotes")
  }
  remotes::install_github("nickbrazeau/goodegg")
  deps <- deps[names(deps) != "goodegg"]
}

# rest of deps
if (any(deps)) {
  install.packages(names(deps)[deps])
}

#......................
# call dependencies
#......................
library(polySimIBD)
library(goodegg)
library(optparse)
library(dplyr)
library(tibble)
library(magrittr)

#++++++++++++++++++++++++++++++++++++++++++
### Functions        ####
#++++++++++++++++++++++++++++++++++++++++++
#............................................................
# utils
#...........................................................
#' @title Internal Function for Extending Migration Matrix in polySimIBD runs
#' @param migmat matrix; original migration matrix
#' @param MigVaryFct integer; amount of time to be spent at home
#' @description Simple liftover function for scaling migration matrix for amount of time that
#' should be spent at home given the original migration matrix (bound, torus, etc.)
#' @returns new migration matrix
#' @export

liftover_migmat <- function(migmat,
                            MigVaryFct) {
  #............................................................
  # checks
  #............................................................
  goodegg::assert_matrix(migmat)
  goodegg::assert_int(MigVaryFct)
  #............................................................
  # core
  #............................................................
  # internal - know this should sum to 100
  # as the MigVaryFct is the percent of time that should be spent at home
  awayperc <- 100 - MigVaryFct
  exchange <- awayperc/sum(migmat)
  migmat[migmat==1] <- exchange
  # now fix diagonal
  diag(migmat) <- MigVaryFct
  #............................................................
  # out
  #............................................................
  return(migmat)
}



#' @title Convert swfsim to IBD
#' @detials Pairwise IBD is from Verity et al 2020. Downsample to allow Ne
#' size to be large but not to overwhelm calculations and to simulate
#' real-life of subset sampling
get_swfsim_2_ibd <- function(swfsim, N, dwnsmplnum){

  # get start and end ind counts for each deme (ie account for when deme size varies)
  # remember, host index is counted as 1:sum(N)
  inds <- lapply(N, function(x){seq(1, x, by = 1)}) # list of inds by deme
  end <- cumsum(sapply(inds, max))
  start <- end + 1 # next start is end + 1, except for last individual
  start <- c(1, start[1:(length(start)-1)])
  # downsample to "N" individuals per deme
  dwnsmpl <- mapply(function(x,y){sample(x:y, size = dwnsmplnum, replace = F)},
                    x = start, y = end, SIMPLIFY = F)
  dwnsmpl <- sort(unlist(dwnsmpl))
  # get combinations
  comb_hosts_df <- t(combn(dwnsmpl, 2))
  comb_hosts_list <- split(comb_hosts_df, 1:nrow(comb_hosts_df))
  # get pairwise IBD
  ibd <- purrr::map_dbl(comb_hosts_list, function(hosts, swf) {
    return(polySimIBD::get_pairwise_bv_ibd(swf = swf, host_index = hosts))
  }, swf = swfsim)
  # tidy up for out
  ret <- tibble::as_tibble(comb_hosts_df, .name_repair = "minimal") %>%
    magrittr::set_colnames(c("smpl1", "smpl2")) %>%
    dplyr::mutate(gendist = as.vector(unlist(ibd)))

  # apply demes
  dms <- as.numeric( cut(dwnsmpl, breaks = c(1,cumsum(N))) ) #inelligent coercion of factor to numeric to represent demes
  demeliftoverx <- tibble::tibble(smpl1 = dwnsmpl,
                                  deme1 = dms)
  demeliftovery <- tibble::tibble(smpl2 = dwnsmpl,
                                  deme2 = dms)
  #......................
  # bring together
  #......................
  ret <- ret %>%
    dplyr::left_join(., demeliftoverx, by = "smpl1") %>%
    dplyr::left_join(., demeliftovery, by = "smpl2")


  return(ret)
}




#............................................................
# main
#...........................................................
#' @title polySimIBD Wrapper
#' @description Based on Maestro Map from 01-simdata
#' @param randseedkey key and random seed for each simulation realization
#' @inheritParams sim_swf nDemes tlim mscale rho pos lambdaCOI - are just the base inputs for sim_swf
#' @returns
run_polySimIBD_swf_2_ibd <- function(randseedkey, nDemes, tlim, rho, pos,
                                     lambdaCOI, mscale, NeVaryMult, baseNe, MigVaryFct, migmat,
                                     dwnsmplnum) {
  #............................................................
  # checks
  #............................................................
  goodegg::assert_int(randseedkey)
  goodegg::assert_int(nDemes)
  goodegg::assert_int(tlim)
  goodegg::assert_int(MigVaryFct)
  goodegg::assert_int(dwnsmplnum)
  goodegg::assert_numeric(lambdaCOI)
  goodegg::assert_numeric(mscale)
  goodegg::assert_numeric(rho)
  goodegg::assert_vector(pos)
  goodegg::assert_vector(NeVaryMult)
  goodegg::assert_matrix(migmat)


  #............................................................
  # setup: liftover migration matrix
  #............................................................
  migmatlifed <- liftover_migmat(migmat = migmat,
                            MigVaryFct = MigVaryFct)

  #............................................................
  # setup: liftover migration matrix
  #............................................................
  Nlifted <- NeVaryMult * baseNe
  Nlifted <- floor(Nlifted)

  #............................................................
  # core
  #............................................................
  # run simulation
  swfsim <- polySimIBD::sim_swf(
    pos = pos,
    rho = rho,
    N = Nlifted,
    m = rep(mscale, nDemes),
    mean_coi = rep(lambdaCOI, nDemes),
    tlim = tlim,
    migr_mat = migmatlifed
  )

  # calculate IBD
  ibddat <- get_swfsim_2_ibd(swfsim = swfsim,
                             N = Nlifted,
                             dwnsmplnum = dwnsmplnum)

  #............................................................
  # Save out to drive and return function complete
  #............................................................
  out <- list(
    swfsim = swfsim,
    ibddat = ibddat
  )
  return(out)
}


#++++++++++++++++++++++++++++++++++++++++++
### Command Line     ####
#++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++
#### parse CL inputs     #####
#++++++++++++++++++++++++++++++++++++++++++
option_list=list(
  make_option(c("-n", "--input"),
              type = "character", default = NULL,
              help = paste("Input filename to read simulation guides from"),
              metavar = "character"),

  make_option(c("-g", "--geodist"),
              type = "character", default = NULL,
              help = paste("Geodist file that correlates to deme locations/distances"),
              metavar = "character"),

  make_option(c("-o", "--output"),
              type = "character", default = NULL,
              help = paste("Output filename to write polySimIBD:sim_swf Simulation Realization and IBD resuls"),
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


#++++++++++++++++++++++++++++++++++++++++++
#### Unpack from CL        #####
#++++++++++++++++++++++++++++++++++++++++++
simguide <- readRDS(opt$input)
randseedkey <- simguide$randseedkey
nDemes <- simguide$nDemes
modname <- simguide$modname
tlim <- simguide$tlim
pos <- simguide$pos[[1]]
rho <- simguide$rho
lambdaCOI <- simguide$lambdaCOI
mscale <- simguide$mscale
NeVaryMult <- simguide$NeVaryMult[[1]]
baseNe <- simguide$baseNe
MigVaryFct <- simguide$MigVaryFct
migmat <- simguide$migmat[[1]]
dwnsmplnum <- simguide$dwnsmplnum
output <- opt$output
geodist <- opt$geodist


### Call Seed for Reproducibility
RNGkind(sample.kind = "Rounding")
set.seed(randseedkey)


#++++++++++++++++++++++++++++++++++++++++++
### Run What you Brung Main ####
#++++++++++++++++++++++++++++++++++++++++++
runout <- run_polySimIBD_swf_2_ibd(randseedkey = randseedkey,
                                   tlim = tlim,
                                   rho = rho,
                                   pos = pos,
                                   nDemes = nDemes,
                                   lambdaCOI = lambdaCOI,
                                   mscale = mscale,
                                   NeVaryMult = NeVaryMult,
                                   baseNe = baseNe,
                                   MigVaryFct = MigVaryFct,
                                   migmat = migmat,
                                   dwnsmplnum = dwnsmplnum)
#......................
# bring in geodist
#......................
geodist <- readRDS(geodist)
runout[[2]] <- runout[[2]] %>%
  dplyr::left_join(., geodist, by = c("deme1", "deme2"))


saveRDS(runout,
        file = output)

# turn warnings back to default
options(warn = defaultwarnings)

