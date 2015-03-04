
# compute the "real" BIC for group-lasso (to check !!)
# ie compute the degree of freedom for each model

modelselection <- function(BICcst, loss.set, beta.refit, beta.avantrefit, size.set, n.nodes, n.samples, mergerefit) { 

  if(n.samples==1 | mergerefit) {
    degreef <- size.set
  } else {
    degreef <- sapply(1:length(beta.refit), FUN=function(ss){
                      tosum <- sapply(1:ncol(beta.refit[[ss]]), FUN=function(pp){
                                      l2norm_refit <- sqrt(sum( beta.refit[[ss]][,pp]^2 ))
                                      if(l2norm_refit > 0){
                                        l2norm_avantrefit <- sqrt(sum( t(beta.avantrefit[[ss]])[,pp]^2 ))
                                        return(1 + (n.samples-1)*l2norm_avantrefit/l2norm_refit)
                                      } else {
                                        return(0) }})
                      return(sum(tosum))}) 
  }
 
  select <- which.min(loss.set + BICcst*degreef*log(n.nodes))  # BICcst useful ? don't know ...

  return(select)

}
