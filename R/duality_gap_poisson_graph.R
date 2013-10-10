duality_gap_poisson_graph <- function(y, paths, beta, graph, param)
{
  ##p <- length(param$loss_weights)
  n <- length(y)
  npaths <- ncol(paths)
  beta0 <- beta
  if(!is.null(npaths) && (npaths > 0)){
    z <- param$loss_weights * (paths %*% beta0)+param$delta
  }else{
    z <- matrix(param$delta, length(y), 1)
    beta0 <- 0
  }
  primal <- sum(z)- sum(y* log(z)) + param$lambda*sum(abs(beta0))
  if(param$lambda==0){
    duality_gap <- Inf
    return()
  }
  ##fprintf('Poisson Loss: %f\n',loss-n*param$delta)
  
  grad1 <- 1.0 - as.matrix(y/z)
  grad2 <- param$loss_weights * grad1
  grad2 <- - grad2/param$lambda
  egp.res <- spams_flipflop.evalPathCoding(grad2, graph, regul=param$regul, resetflow=TRUE, eval_dual=TRUE,
                                  pos=param$pos)
  dual_norm <- egp.res$dual

  grad1 <- grad1/max(dual_norm,1)

  thrs <- 1.0 - y/param$delta
  fenchel <- 0
  for(ii in 1:n){
    if (grad1[ii] <= thrs[ii]){
      fenchel <- fenchel + (-param$delta + y[ii]*log(param$delta))
    }else{
      if (grad1[ii] <= 1.0){
        fenchel <- fenchel + (-param$delta*grad1[ii] - y[ii] + y[ii]*log(y[ii]/(1e-15+1.0-grad1[ii])))
      }else{
        fenchel <- fenchel + Inf
      }
    }
  }

  dual <- -fenchel
  duality_gap <- (primal-dual)/abs(primal)
  return(list(duality_gap=duality_gap, primal=primal, dual=dual))
}
