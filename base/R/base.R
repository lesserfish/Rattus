
# system is the system dynamics. It is a function that receives the following inputs:
# function(t, x, params).
# It output the derivative of x at time t under the specified parameters

# perturbation is the change in the system parameters. It receives a list as input with the current context.
# The context list contains an element 'iteration' consisting of the current iteration level.


apply_cperturbation <- function(system, x0, perturbation, end_condition, params = list(), t0=0, tf=1, teps=0.1)
{
  it <- 1
  context <- list('params' = params, 'iteration' = it)
  trajectories <- list()

  while(!end_condition(context)) {
    # Calculate Trajectory

    times <- seq(t0, tf, teps)

    current_trajectory <- deSolve::ode(func = system, y = x0, times=times, parms=params)

    tmedian <- apply(current_trajectory[,-1], MARGIN=2, median)
    xf <- current_trajectory[dim(current_trajectory)[1], -1]
    meanf <- mean(xf)

    trajectories[[it]] <- list('trajectory' = current_trajectory,
                               'final_point' = xf,
                               'meanf' = meanf,
                               'median' = tmedian,
                               'context' = context)

    context <- list('params' = params,
                    'iteration' = it ,
                    'trajectories' = trajectories)

    # Update parameters
    params <- perturbation(context)
    it <- it + 1
  }

  trajectories[['total_iterations']] <- it

  return(trajectories)
}
