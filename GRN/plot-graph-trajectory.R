# Load the igraph package
library(igraph)
library(RColorBrewer)
library(deSolve)
library(animation)

set.seed(123)

color <- function(x, N = 100) {
  my_palette <- colorRampPalette(c("black", "yellow", "orange", "red", "white"))
  output <- my_palette(N+1)[floor(N * x)+1]
  return(output)
}

Michaelis_Menten <- function(t, x, params)
{
  A <- params$A
  f <- params$f
  h <- params$h
  B <- params$B
  outer <- params$outer
  dx <-  B*(x * outer)^f + (x^h / (x^h + 1)) %*% A
  #dx <-  B*(x)^f + (x^h / (x^h + 1)) %*% A
  return(list(dx))
}

find_outer_edges <- function(A)
{
  output <- rep(1, dim(A)[1])
  incoming_edges <- colSums(A)
  output[incoming_edges == 0] <- 0
  return(output)
}

load("~/Documents/Work/SysBio2023/CHIKV Data/scripts/ko04657.RData")
load("~/Documents/Work/SysBio2023/Presentations/networks/ko04062.RData")
A <- as.matrix(igraph::as_adjacency_matrix(full_graph, names=FALSE))

#set.seed(1231232)
#A <- as.matrix(igraph::as_adjacency_matrix(igraph::sample_gnp(10, 0.1, directed=TRUE)))

outer <- find_outer_edges(A)

graph <- igraph::graph_from_adjacency_matrix(A)
V(graph)$color <- color(1.0*(1 - outer))
V(graph)$border.color <- "white"
layout <- igraph::layout.auto(graph)

par(bg="black")

plot(graph, layout=layout, vertex.frame.color="white", vertex.size=7, edge.arrow.size=0.4, vertex.label.color = "white", vertex.label.cex = 0.5)


ts <- deSolve::ode(func = Michaelis_Menten,
                   y = 2*(1 - outer),
                   #y = rep(2, dim(A)[1]),
                   times = seq(1, 50),
                   parms = list("A" = A, "f" = 1, "h" = 1, 'B' = -1.0, 'outer'=outer))

maxval <- 3.0

dts <- diff(ts[,-1])
rge <- which(rowSums(abs(dts)) < 0.05)[1] + 5

saveGIF({
  for(i in 1:rge) {
    par(bg="black")
    g <- graph
    values <- ts[i, -1]
    values <- pmax(0, pmin(values, maxval))
    values <- (1/maxval) * values

    V(g)$color <- color(values)
    plot(g, layout = layout, vertex.size=7, edge.arrow.size=0.4, vertex.frame.color="white", vertex.label.color = "black", vertex.label.cex = 0.5)

  }
}, movie.name = "~/Documents/Work/SysBio2023/Presentations/network-flow-tmp.gif", interval = 0.1)

write.csv(ts, "/tmp/a.csv")

