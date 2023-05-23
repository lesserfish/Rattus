# Try with real Network
library(igraph)


par(bg="black")
plot(NA, NA, xlim=c(0, 1), ylim = c(0, 6))
axis(1, col="white", col.ticks="white", col.axis="white", cex.axis=1)
axis(2, col="white", col.ticks="white", col.axis="white", cex.axis=0.8)


setwd('~/Documents/Code/Rattus/GRN/')
load("./IL17.RData")
G <- full_graph
A <- as.matrix(igraph::as_adjacency_matrix(G, names=FALSE))
n <- dim(A)[1]

itcount <- 1

for(iii in 1:itcount) {

output <- matrix(NA, nrow=itcount, ncol=n)

A <- as.matrix(igraph::as_adjacency_matrix(G, names=FALSE))
A[A != 0] <- 2

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
  tf = 400
  teps = 0.68


  result <- apply_cperturbation(system = Michaelis_Menten,
                                 x0 = x0,
                                 perturbation = perturbation,
                                 end_condition = end_condition,
                                 params = parameters,
                                 t0 = t0,
                                 tf= tf,
                                 teps=teps,
                                 log=FALSE)


for(i in 1:n)
{
  gene_trajectory <- result[,i] # Vector of end states of gene i
  lines(seq(0, 1, length.out = length(gene_trajectory)), y = gene_trajectory, col = rgb(1, 0, 0, alpha = 0.5))
}

# Plot overall
  medians <- apply(result, MARGIN = 1, median)
  lines(seq(0, 1, length.out = length(medians)), y = medians, col = rgb(1, 1, 1, alpha = 1.0), lwd=2)

means <- apply(result, MARGIN = 1, mean)
lines(seq(0, 1, length.out = length(means)), y = means, col = rgb(0, 0, 1, alpha = 1.0), lwd=2)


}

