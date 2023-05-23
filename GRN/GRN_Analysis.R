library(igraph)

# Generate random graph to use as a demo
n <- 100
G <- igraph::sample_pa(n)
A <- 2*as.matrix(igraph::as_adjacency_matrix(G))

#Plot G
{
  V(G)$size <- 5
  V(G)$color <- "#722f37"
  V(G)$label.cex <- 0.5
  E(G)$arrow.mode <- 0
  plot(G, size=0)
}

Michaelis_Menten <- function(t, x, params)
{
  A <- params$A
  f <- params$f
  h <- params$h

  dx <- -x^f + A %*% (x^h / (x^h + 1))

  return(list(dx))
}

end_condition <- function(context)
{
  it <- context$iteration
  # Stop the perturbation once we have removed all nodes
  if(it == n){
    return(TRUE)
  }
  return(FALSE)
}


removal_order <- sample(1:n, replace=FALSE)

perturbation <- function(context)
{
  params <- context$params
  it <- context$iteration

  A <- params$A

  # Select which node will be removed
  node <- removal_order[it]
  # Remove outgoing edges
  A[node,] <- 0

  # Remove incoming edges
  A[, node] <- 0

  params[['A']] <- A
  return(params)
}


parameters <- list('f' = 1, 'h'=2, 'A'=A)
x0 = rep(2, n)
t0 = 0
tf = 100
teps = 1


result <- apply_cperturbation2(system = Michaelis_Menten,
                              x0 = x0,
                              perturbation = perturbation,
                              end_condition = end_condition,
                              params = parameters,
                              t0 = t0,
                              tf= tf,
                              teps=teps,
                              log=FALSE)

# Get results
plot(x=1:length(result), y = result, type='l', col='red')




