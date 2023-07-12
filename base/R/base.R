
# system is the system dynamics. It is a function that receives the following inputs:
# function(t, x, params).
# It output the derivative of x at time t under the specified parameters

# perturbation is the change in the system parameters. It receives a list as input with the current context.
# The context list contains an element 'iteration' consisting of the current iteration level.


apply_cperturbation <- function(system, x0, perturbation, end_condition, params = list(), t0=0, tf=1, teps=0.1, log=FALSE)
{
  it <- 1
  context <- list('params' = params, 'iteration' = it)

  end_state <- matrix(nrow = 0, ncol=length(x0))

  while(!end_condition(context)) {
    # Calculate Trajectory

    if(log) {
      print(paste("Current iteration: ", it, sep=""))
    }
    times <- seq(t0, tf, teps)

    current_trajectory <- deSolve::ode(func = system, y = x0, times=times, parms=params)

    xf <- current_trajectory[dim(current_trajectory)[1], -1]
    end_state <- rbind(end_state, xf)

    context <- list('params' = params,
                    'iteration' = it )

    # Update parameters
    params <- perturbation(context)
    it <- it + 1
  }

  return(end_state)
}

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

node_removal_end_condition <- function(context)
{
  it <- context$iteration
  # Stop the perturbation once we have removed all nodes
  if(it == edge_count){
    return(TRUE)
  }
  return(FALSE)
}

node_removal_perturbation <- function(A)
{
  n <- dim(A)[1]
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
  return(perturbation)
}

node_removal_end_condition <- function(context)
{
  n <- dim(context$params$A)[1]
  it <- context$iteration
  # Stop the perturbation once we have removed all nodes
  if(it == n){
    return(TRUE)
  }
  return(FALSE)
}

node_removal_parameters = function(A)
{
  outer <- find_outer_edges(A)
  parameters <- list('f' = 1, 'h'=2, 'A'=A, 'B'=-1, 'outer'=outer)
  return(parameters)
}

perturbation_analysis <- function(M,
                      x0,
                      t0=0,
                      tf=40,
                      teps = 0.1,
                      perturbation_function = node_removal_perturbation,
                      end_condition = node_removal_end_condition,
                      parameter_function = node_removal_parameters,
                      iterations = 1,
                      dynamics = Michaelis_Menten)
{
  n <- dim(M)[1]
  output <- matrix(NA, nrow=iterations, ncol=n)
  for(iteration in 1:iterations) {

    A <- M


    this_perturbation <- perturbation_function(A)
    parameters <- parameter_function(M)

    result <- apply_cperturbation2(system = dynamics,
                                   x0 = x0,
                                   perturbation = this_perturbation,
                                   end_condition = end_condition,
                                   params = parameters,
                                   t0 = t0,
                                   tf= tf,
                                   teps=teps,
                                   log=FALSE)

    output[iteration, ] <- result
  }
  return(output)
}

plot_resilience_drop <- function(data)
{
  par(bg="black")
  max_y <- data[1, 1]
  plot(NA, NA, xlim=c(0, 1), ylim=c(0,  max_y))
  axis(1, col="white", col.ticks="white", col.axis="white", cex.axis=1)
  axis(2, col="white", col.ticks="white", col.axis="white", cex.axis=0.8)

  iterations <- dim(data)[1]
  for(iteration in 1:iterations)
  {
    result <- data[iteration, ]
    lines(seq(0, 1, length.out = length(result)), y = result, col = rgb(1, 0, 0, alpha = 0.5))
  }

  outputmed <- apply(data, MARGIN = 2, FUN = median)
  lines(seq(0, 1, length.out = length(outputmed)), y = outputmed, col = rgb(1, 1, 1, alpha = 1.0), lwd=2)
  legend("topright", legend = c("Gene Expression", "Median"), col = c("red", "white"), lty = 1, fill = c("red", "white"), text.col = "white")
}
