#..............................................................
# Purpose of this script is to find optimal
# lambda values that give us a mean COI of 1, 2, 3,4 for
# our models
#..............................................................
set.seed(48)
minLambda <- function(x, y){
  yhat <- x / (1 - exp(-x))
  return( (y-yhat)^2 )
}

lamdba_optimdf <- tibble::tibble(
  meancoi = c(1, 2, 3, 4),
  start = c(1e-5, 2.3, 3.1, 4.2)
)

lamdba_optimdf$optimpar <- purrr::pmap(lamdba_optimdf, function(meancoi, start){
  ret <- optim(par = start, fn = minLambda, method = "BFGS", y = meancoi)
  return(ret$par)
})

# make sure positive for first optim
# it is all approximately the same (as it approximating 0, L'hospital)

# save out
dir.create("01-simdata/sim_params/", recursive = T)
saveRDS(object = abs(unlist(lamdba_optimdf$optimpar)),
        file = "01-simdata/sim_params/optim_lambda.RDS")
