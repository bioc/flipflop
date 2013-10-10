max_lambda <- function(y, graph, delta, loss_weights, pos)
{
   ## X <- as(diag(loss_weights), 'CsparseMatrix')
   grad1  <- 1.0 - y/delta
   grad2 <-  loss_weights * grad1
   grad2  <-  -grad2
   regul <- 'graph-path-conv2'
   egp.res <- spams_flipflop.evalPathCoding(grad2, graph, regul=regul, 
                                   resetflow=TRUE, eval_dual=TRUE,
                                   pos=pos)
   return(egp.res)
}
