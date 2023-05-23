# Try with real Network
library(igraph)


par(bg="black")
plot(NA, NA, xlim=c(0, 1), ylim = c(0, 20))
axis(1, col="white", col.ticks="white", col.axis="white", cex.axis=1)
axis(2, col="white", col.ticks="white", col.axis="white", cex.axis=0.8)


setwd('~/Documents/Code/Rattus/GRN/')
load("./IL17.RData")
#load("~/Documents/Work/SysBio2023/Presentations/networks/ko04062.RData")
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
  B <- params$B
  outer <- params$outer
  dx <-  B*(x * outer)^f + (x^h / (x^h + 1)) %*% A
  return(list(dx))
}

find_outer_edges <- function(A)
{
  output <- rep(1, dim(A)[1])
  incoming_edges <- colSums(A)
  output[incoming_edges == 0] <- 0
  return(output)
}

outer <- find_outer_edges(A)

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
    outer <- params$outer

    # Select which node will be removed
    node <- removal_order[it]
    # Remove outgoing edges
    A[node,] <- 0
    outer[node] <- 1 # If we remove an outer node, it's expression can (and will) drop
    # Remove incoming edges
    A[, node] <- 0

    params[['A']] <- A
    params[['outer']] <- outer
    return(params)
  }


  parameters <- list('f' = 1, 'h'=2, 'A'=A, 'B'=-1, 'outer'=outer)
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
  lines(seq(0, 1, length.out = length(gene_trajectory)), y = gene_trajectory, col = "#6b1c23", lwd=0.7)
}

# Plot overall
#  medians <- apply(result, MARGIN = 1, median)
#  lines(seq(0, 1, length.out = length(medians)), y = medians, col = rgb(1, 1, 1, alpha = 1.0), lwd=2)

means <- apply(result, MARGIN = 1, mean)
lines(seq(0, 1, length.out = length(means)), y = means, col = rgb(1, 0, 0, alpha = 1.0), lwd=1.5)

legend("topright", legend = c("Gene expression", "Mean"), col = c("#6b1c23", "red"), lty = 1, fill = c("#6b1c23", "red"), text.col = "white")
title(main = "IL-17 Signaling path drop rate", col.main="white")

}

