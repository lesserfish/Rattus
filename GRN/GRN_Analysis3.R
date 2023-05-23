# Try with real Network

library(R.matlab)
library(Matrix)

data <- readMat("~/Documents/Work/SysBio2023/NuRsE/simulation_real/R_simulation_code/R_real_data/TECB.mat")
A <- data$A
n <- dim(A)[1]
G <- igraph::graph_from_adjacency_matrix(A)
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
  if(it == 40){
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
tf = 10
teps = 1

out <- deSolve::ode(y = rep(2, n), times = seq(0, 100, 1), func = Michaelis_Menten, parms = parameters)
mean(out[101, -1])

result <- apply_cperturbation2(system = Michaelis_Menten,
                              x0 = x0,
                              perturbation = perturbation,
                              end_condition = end_condition,
                              params = parameters,
                              t0 = t0,
                              tf= tf,
                              teps=teps,
                              log=TRUE)

# Get results
plot(x=1:length(result), y = result, type='l', col='red')



