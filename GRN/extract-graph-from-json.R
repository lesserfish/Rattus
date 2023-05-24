library(jsonlite)
library(igraph)
library(R.matlab)

directory <- "/home/lesserfish/Documents/Work/SysBio2023/Presentations/networks" # A pasta onde o arquivo esta
jsonfile <- "ko05323.json" # O nome do arquivo Json
name <- "true name" # Escrever algum codigo para transformar de codigo a nome. Ex: ko04657 -> IL17 Signaling Path. Possivelmente usar o arquivo pathways.txt


filepath <- paste(directory, "/", jsonfile, sep="")
data <- fromJSON(filepath)

nodes <- data$nodes
edges <- data$edges

G <- igraph::make_empty_graph()

node_dict <- list()

for(i in 1:length(nodes[[2]]))
{
  id <- nodes[[1]][i]
  node <- nodes[[2]][i]
  node_dict[[id]] <- node

  G <- igraph::add.vertices(G, 1, name=node)
}

for(i in 1:length(edges[[2]]))
{
  sourceid <- edges[[1]][i]
  targetid <- edges[[2]][i]
  type <- edges[[3]][i]
  subtype <- edges[[4]][i]

  source <- node_dict[[sourceid]]
  target <- node_dict[[targetid]]

  G <- igraph::add.edges(G, c(source, target), subtype=subtype)

}

G <- igraph::delete.vertices(G, which(igraph::degree(G) == 0))


V(G)$size <- 5
V(G)$color <- "#722f37"
V(G)$label.cex <- 0.5
E(G)$arrow.mode <- 0
plot(G, size=0)  # Uncomment this if you want to plot the graph


rdataname <- sub(".json", ".RData", jsonfile)
rdatafilepath <- paste(directory, "/", rdataname, sep="")
save(G, name, file = rdatafilepath)

