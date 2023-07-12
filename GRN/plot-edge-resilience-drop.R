# Try with real Network
library(igraph)


par(bg="black")
plot(NA, NA, xlim=c(0, 1), ylim=c(0, 5))
axis(1, col="white", col.ticks="white", col.axis="white", cex.axis=1)
axis(2, col="white", col.ticks="white", col.axis="white", cex.axis=0.8)


setwd('~/Documents/Code/Rattus/GRN/')
load("./IL17.RData")
G <- full_graph
A <- as.matrix(igraph::as_adjacency_matrix(G, names=FALSE))
edge_count <- length(which(as.vector(A) != 0))
node_count <- dim(A)[[1]]
itcount <- 200

output <- matrix(NA, nrow=itcount, ncol=edge_count)

for(iii in 1:itcount) {
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
    if(it == edge_count){
      return(TRUE)
    }
    return(FALSE)
  }

  vA <- as.vector(A)
  edges <- which(vA != 0)
  edges <- edges[sample(1:length(edges), replace=FALSE)]
  perturbation <- function(context)
  {
    params <- context$params
    it <- context$iteration

    A <- params$A
    outer <- params$outer

    # Select which node will be removed
    vA <- as.vector(A)
    e <- edges[it]
    vA[e] <- 0
    A <- matrix(vA, nrow=nrow(A))

    # TODO: Check if any outer node is now empty

    params[['A']] <- A
    params[['outer']] <- outer
    return(params)
  }


  parameters <- list('f' = 1, 'h'=2, 'A'=A, 'B'=-1, 'outer'=outer)
  x0 = rep(2, node_count)
  t0 = 0
  tf = 400
  teps = 0.68


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
  lines(seq(0, 1, length.out = length(result)), y = result, col = rgb(0, 1.0, 0, alpha = 0.5))

  output[iii, ] <- result
}

outputm <- apply(output, MARGIN = 2, FUN = mean)
lines(seq(0, 1, length.out = length(outputm)), y = outputm, col = rgb(0, 0.0, 1.0, alpha = 1.0), lwd=2)

outputmed <- apply(output, MARGIN = 2, FUN = median)
lines(seq(0, 1, length.out = length(outputmed)), y = outputmed, col = rgb(1, 1, 1, alpha = 1.0), lwd=2)

legend("topright", legend = c("Mean", "Median"), col = c("blue", "white"), lty = 1, fill = c("blue", "white"), text.col = "white")

title(main = "IL-17 Signaling path edge drop rate", col.main="white")
