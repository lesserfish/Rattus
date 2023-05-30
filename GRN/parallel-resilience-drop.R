library(igraph)
library(foreach)
library(doParallel)


NUM_CORES <- 64   # Number of cores
ITER = 1000      # Number of iterations per network
TEXT_COLOR = "white"
BACKGROUND_COLOR = "black"
GENE_COLOR = rgb(1, 0, 0, alpha = 0.5)
MEDIAN_COLOR = rgb(1, 1, 1, alpha = 1.0)
basedir = "~/Documents/Work/SysBio2023/Kegg/Kegg_Rdata/Rdata/"

# Get the list of available networks

network_fp = paste(basedir, "/networks", sep="")
file <- file(network_fp, "r")
networks <- readLines(file)
close(file)

# How many networks we have
COUNT <- length(networks)

# Perturbation Analysis functions
apply_cperturbation2 <- function(system, x0, perturbation, end_condition, params = list(), t0=0, tf=1, teps=0.1, log=FALSE)
{
  it <- 1
  context <- list('params' = params, 'iteration' = it)
  means <- c()
  
  while(!end_condition(context)) {
    # Calculate Trajectory
    
    if(log) {
      print(paste("Current iteration: ", it, sep=""))
    }
    times <- seq(t0, tf, teps)
    
    current_trajectory <- deSolve::ode(func = system, y = x0, times=times, parms=params)
    
    xf <- current_trajectory[dim(current_trajectory)[1], -1]
    meanf <- mean(xf)
    
    means <- c(means, meanf)
    
    context <- list('params' = params,
                    'iteration' = it )
    
    # Update parameters
    params <- perturbation(context)
    it <- it + 1
  }
  
  return(means)
}

# The Model we are using
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

find_outer_nodes <- function(A)
{
  output <- rep(1, dim(A)[1])
  incoming_edges <- colSums(A)
  output[incoming_edges == 0] <- 0
  return(output)
}

Analyze_Network <- function(network_id)
{
  # Loads the network onto memory

  network_code <- networks[network_id]
  network_dir <- paste(basedir, "/", network_code, sep="")
  network_data_fp <- paste(network_dir, "/", network_code, ".RData", sep="")
  load(network_data_fp)
  network_name = name[[1]]
  
  #
  message <- sprintf(": Began analysis of network %s (%s)", network_code, network_name)
  cat(message, "\n")
  
  # Convert graph to Matrix
  A <- as.matrix(igraph::as_adjacency_matrix(G, names=FALSE))
  n <- dim(A)[1]
  
  output <- matrix(NA, nrow=ITER, ncol=n)
  
  for(iteration in 1:ITER){
    
    # Remove any perturbation that was already caused to the matrix
    A <- as.matrix(igraph::as_adjacency_matrix(G, names=FALSE))
    A[A != 0] <- 2
    
    # Calculate the outer nodes
    outer <- find_outer_nodes(A)
    
    # Specify an end condition to the perturbations
    
    end_condition <- function(context)
    {
      it <- context$iteration
      # Stop the perturbation once we have removed all nodes
      if(it == n){
        return(TRUE)
      }
      return(FALSE)
    }
    # Establish the order in which the nodes will be removed
    removal_order <- sample(1:n, replace=FALSE)
    
    # Function that perturbates the network
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
    # Set the parameters to call the cperturbation_analysis
    parameters <- list('f' = 1, 'h'=2, 'A'=A, 'B'=-1, 'outer'=outer)
    x0 = rep(2, n)
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
    # Store the result to output
    output[iteration, ] <- result
  }
  # Calculate the median of the graphs
  median <- apply(output, MARGIN = 2, FUN = median)
  
  # Calculate Median integral
  base <- (1/length(output))
  integral <- sum(base * output)
  
  # Plot results
  plot_fp <- paste(network_dir, "/", network_code, "_node_resilience.png", sep="")
  png(plot_fp)
  
  ymax <- output[1, 1]
  
  par(bg=BACKGROUND_COLOR, col=TEXT_COLOR)
  plot(NA, NA, 
       xlim=c(0, 1), 
       ylim=c(0, ymax), 
       xlab="Proportion of removed genes", 
       ylab="Gene expression", 
       col.axis=TEXT_COLOR, 
       col=TEXT_COLOR, 
       col.lab=TEXT_COLOR, 
       col.main=TEXT_COLOR, 
       col.sub=TEXT_COLOR, 
       col.axis=TEXT_COLOR)
  
  for(iteration in 1:ITER){
    result <- output[iteration, ]
    lines(seq(0, 1, length.out = length(result)), y = result, col = GENE_COLOR)
  }
  lines(seq(0, 1, length.out = length(median)), y = median, col = MEDIAN_COLOR, lwd=2)
  legend("topright", legend = c("Gene Expression", "Gene Expression Median"), col = c(GENE_COLOR, MEDIAN_COLOR), lty = 1, fill = c(GENE_COLOR, MEDIAN_COLOR), text.col = TEXT_COLOR)
  title(main=network_name, col.main=TEXT_COLOR)
  dev.off()

  # Normalize results
  
  noutput <- output/output[1]
  nmedian <- apply(noutput, MARGIN = 2, FUN = median)
  
  nbase <- (1/length(noutput))
  nintegral <- sum(base * noutput)
  
  # Plot results
  plot_fp <- paste(network_dir, "/", network_code, "_node_resilience_normalized.png", sep="")
  png(plot_fp)
  
  ymax <- 1
  
  par(bg=BACKGROUND_COLOR, col=TEXT_COLOR)
  plot(NA, NA, 
       xlim=c(0, 1), 
       ylim=c(0, ymax), 
       xlab="Proportion of removed genes", 
       ylab="Gene expression", 
       col.axis=TEXT_COLOR, 
       col=TEXT_COLOR, 
       col.lab=TEXT_COLOR, 
       col.main=TEXT_COLOR, 
       col.sub=TEXT_COLOR, 
       col.axis=TEXT_COLOR)
  
  for(iteration in 1:ITER){
    result <- noutput[iteration, ]
    lines(seq(0, 1, length.out = length(result)), y = result, col = GENE_COLOR)
  }
  lines(seq(0, 1, length.out = length(nmedian)), y = nmedian, col = MEDIAN_COLOR, lwd=2)
  legend("topright", legend = c("Gene Expression", "Gene Expression Median"), col = c(GENE_COLOR, MEDIAN_COLOR), lty = 1, fill = c(GENE_COLOR, MEDIAN_COLOR), text.col = TEXT_COLOR)
  title(main=network_name, col.main=TEXT_COLOR)
  dev.off()
  
  # Save results to disk
  output_fp <- paste(network_dir, "/", network_code, "_gene_analysis.RData", sep="")
  save(G, name, output, median, integral, nmedian, nintegral, file = output_fp)
  message <- sprintf(": Finished analysis of network %s (%s)", network_code, network_name)
  cat(message, "\n")
  
}

# Perform analysis in multiple cores

# Register cores
cl <- makeCluster(NUM_CORES)
registerDoParallel(cl)

result <- foreach(i = 1:COUNT) %dopar% {
  Analyze_Network(i)
}

stopCluster(cl)

