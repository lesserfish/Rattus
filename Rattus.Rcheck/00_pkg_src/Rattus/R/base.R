
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

