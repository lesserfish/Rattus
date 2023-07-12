library(igraph)
load("~/Documents/Work/SysBio2023/KEGG_DATABASE_FULL/ko00360/ko00360.RData")
M <- as.matrix(igraph::as_adjacency_matrix(G, names=FALSE))
n <- dim(M)[1]
x0 = rep(2, n)

data <- perturbation_analysis(M = M, x0 = x0, iterations = 10)

plot_resilience_drop(data)
