#---------------------------------#
#---------- GROUP LASSO ----------#
#----------------l1/2-------------#

# *** refit transcript abundances under poisson loss without penalization *** given known design (transcripts)
# INPUT:        - count vector y per samples : nxT
#               - design X, ie. found transcript structures : nxp
#               - initialization beta0, here the shrinked values of abundances : pxT
#               - delta smoothing loss
#               - optimization parameters iterpoisson and tolpoisson
#               - 'solver_refit type of solver'
# OUTPUT:       - beta.refit the values of abundances for given regularization

grouplassoPoisson <- function(y, X, beta0, delta, iterpoisson, tolpoisson, solver_refit, verbosepath, n.samples, lambda){

#   if(solver_refit=='OLD'){
#      beta.refit <- spams_flipflop_flipflop.solverPoisson(y=y, X=X, beta0=beta0, weights=rep(0, length(beta0)), delta=delta, max_iter=iterpoisson, tol=tolpoisson)
#   }
 
   W0 <- beta0
   if(solver_refit=='NEW'){
      seqdelt <- 10^seq(0, log10(delta), length=(-floor(log10(delta))+1))
      seqdelt <- c(seqdelt, delta)  ## important OR NOT ??????
      for(deltaloop in seqdelt){
         beta.refit <- spams_flipflop.fistaFlat(Y=y, X=X, W0=W0, 
                                       loss='poisson', regul='l1l2', delta=deltaloop, pos=TRUE, 
                                       tol=tolpoisson, max_it=iterpoisson, ista=TRUE, linesearch_mode=2, L0=1e-5, verbose=verbosepath,
                                       size_group=n.samples, lambda1=lambda)
         W0 <- beta.refit
      }
   }

   return(beta.refit)

}
