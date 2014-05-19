
duality_gap_poisson <- function(y, X, w, delta){
  
  z <- X %*% w + delta
  n <- nrow(z)

  loss <- sum(z) - sum(y * log(z))
  ##print(sprintf('Poisson Loss: %f', loss - n * delta))

  vec1 <- colSums(X)
  tmp <- y/z
  vec2 <- t(X) %*% tmp
  vec <- vec1 / vec2
  gamma <- min(min(vec), 1.0)
  nu <- rep(1, n) - gamma * tmp
  dual <- sum(y) - sum(y * log( z/gamma)) + delta * sum(nu)
  duality_gap <- (loss-dual) / abs(loss)
  ##print(sprintf('Relative duality gap: %f', duality_gap))
  return(list(dg=duality_gap, loss=loss))
}

