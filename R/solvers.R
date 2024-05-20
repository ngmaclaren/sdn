#' @import graphics stats parallel deSolve
library(graphics)
library(stats)
library(parallel)
library(deSolve)

#' Simulate ODEs
#' 
#' @description Simulate a system of ODEs across a range of a control parameter.
#'
#' @param rng The range of values the control parameter takes. 
#' @param varname The name of the control parameter. Must match a parameter in params list.
#' @param func The dynamical model, a function returning a derivative
#' @param initialvalue A vector or scalar of initial value(s).
#' @param params A list of function parameters, including the adjacency matrix if necessary.
#' @param control A list including some subset of nsamples, spacing, deltaT, times, ncores, or silent, depending on the application.
#' @param kind To solve with or without noise
#' @param allsamples Not currently used
#' 
#' @details Solves an ODE or SDE across a range of parameter values. If `kind = "ode"`, requires deSolve. 
#'
#' @return An L x N matrix of final values, where L is the length of `rng` and N is the length of `initialvalue`.
#' @examples
#' library(graphics)
#' library(stats)
#' library(deSolve)
#' library(localsolver)
#' 
#' ## Generate an adjacency matrix
#' library(igraph)
#' g <- largest_component(sample_gnm(10, 20, directed = FALSE, loops = FALSE))
#' A <- as_adj(g, "both", sparse = FALSE)
#' N <- vcount(g)
#' 
#' ## solve_in_range() assumes mclapply() is available, so may not work on Windows
#' library(parallel)
#' ncores <- detectCores()-1
#' 
#' params <- c(.doublewell, list(A = A))
#' params$sigma <- 0.1
#' 
#' control <- list(times = 0:10, deltaT = 0.01, ncores = ncores)
#' simtimes <- seq(control$times[1], control$times[length(control$times)], by = control$deltaT) # for comparison between ode() and sde()
#' 
#' X.ode <- solve_in_range(params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode")
#' X.sde <- solve_in_range(params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "sde")
#' 
#' xlim <- range(params$Ds)
#' ylim <- range(c(X.ode, X.sde))
#' 
#' matplot(params$Ds, X.ode, xlim = xlim, ylim = ylim, xlab = "D", ylab = "x", lty = 1, col = 1, lwd = 1, type = "l")
#' matlines(params$Ds, X.sde, lty = 1, col = 2, lwd = 1)
#' @export
solve_in_range <- function(
                           rng, # the range of the parameter value, given as a vector
                           varname, # the name of the control parameter in the model, a string
                           func, # the dynamical model, returning a derivative with or without noise
                           initialvalue, # pass as vec or scalar... needs to be done outside
                           params = list(), # for the func
                           control = list(), # nsamples, spacing, deltaT, times, ncores, silent
                           kind = "ode", # or "sde"
                           allsamples = FALSE # returns 3D array
                           ) {
                                        # If mclapply won't run, try with ncores = 1
    if(!("ncores" %in% names(control))) control$ncores <- 1
    if(!("silent" %in% names(control))) control$silent <- FALSE # this doesn't seem to do much
    if(!("times" %in% names(control))) control$times <- seq(0, 10, by = 1)

    if(kind == "ode") {
        results <- mclapply(
            rng, function(val) {
                params[[varname]] <- val
                ode(initialvalue, control$times, func, params)
            }, mc.cores = control$ncores, mc.silent = control$silent
        )
                                        # Treat the final value of x_i as x_i^*
        result <- do.call(rbind, lapply(results, function(df) df[nrow(df), -1]))
        
    } else if(kind == "sde") {
        results <- mclapply(
            rng, function(val) {
                params[[varname]] <- val 
                sde(initialvalue, control$times, func, params, control)
            }, mc.cores = control$ncores, mc.silent = control$silent
        )
                                        # Treat the final value of x_i as x_i^*
        result <- do.call(rbind, lapply(results, function(df) df[nrow(df), -1]))
    }
    rownames(result) <- seq_len(nrow(result))
    return(result)
}

### SDEs

#' Determine the number of time steps for a simulation
#'
#' Uses the desired time span and the simulation time step to compute how many time steps are required for the simulation.
#' @param times a sequence of time steps, in user units
#' @param maxT The last time step (in user units)
#' @param minT The first time step (in user units)
#' @param deltaT The size of the timestep for simulation
#' @return A scalar
#' @examples
#' determine_ntimesteps(times = 0:10, deltaT = 0.01)
#' @export
determine_ntimesteps <- function(times = numeric(), maxT = 1, minT = 0, deltaT = 0.01) {
                                        # Requires times be at least length 2.
                                        # Not sure that's what I want.
    if(length(times > 1)) {
        minT <- times[1]
        maxT <- times[length(times)]
    }
    
    (maxT - minT)/deltaT
}
    
#' Preallocate noise for sde()
#'
#' Generates a ntimesteps x nnodes matrix of Gaussian random numbers to use for solving SDEs.
#'
#' @param sigma The standard deviation of the noise process
#' @param ntimesteps The number of simulation time steps (probably the output of determine_ntimesteps())
#' @param nnodes The number of SDEs (i.e., number of nodes in a network)
#' @return A matrix
#' @examples
#' ntimesteps <- determine_ntimesteps(times = 0:10, deltaT = 0.01)
#' W <- preallocate_noise(0.01, ntimesteps, 5)
#' @export
preallocate_noise <- function(sigma, ntimesteps, nnodes) {
                                        # hard-coding Gaussian noise
                                        # Would like to allow for Levy flights, eventually (power-law noise?)
                                        # See here: https://en.wikipedia.org/wiki/L%C3%A9vy_process
                                        # And here: https://en.wikipedia.org/wiki/L%C3%A9vy_flight
                                        # Looks like options are Wiener/Gaussian (below), Poisson, and Cauchy.
                                        # It looks like determining the matrix W would be more complicated for either a Poisson or Cauchy process, though potentially worth doing.
    if(length(sigma) == 1) {
        matrix(
            rnorm(nnodes*ntimesteps, mean = 0, sd = sigma), ncol = nnodes, nrow = ntimesteps
        )
    } else {
        stopifnot(length(sigma) == nnodes) # require σ_i ∀ i, but not necessarily unique
        do.call(
            cbind,
            lapply(sigma, function(s) rnorm(ntimesteps, mean = 0, sd = s))
        )
    }
}

#' Simulate stochastic differential equations on networks
#'
#' Simulate stochastic differential equations using the Euler-Maruyama method. Output mimics the output of deSolve's ode().
#' @param initialvalue The inital value of each variable/node, scalar or vector
#' @param times A sequence, at least the first and last time in user units, but e.g. 0:10 works and is easy.
#' @param func The dynamical model, in deSolve's format. Should return a deterministic derivative. 
#' @param parms A list. Using deSolve's naming convention, the model parameters (including the adjacency matrix, if using). Must include an element called "sigma" for the standard deviation of the noise process.
#' @param control A list. Must include an element called "deltaT".
#' @return A data frame
#' @details Returns a data frame that looks like deSolve's ode() output: a column of time steps, then a column of values at each timestep for each variable in the model.
#' @examples
#' library(graphics)
#' library(stats)
#' library(deSolve)
#' library(localsolver)
#' 
#' .growth <- list(xinit = 0.1, r = 1, K = 5, sigma = 0.1) # Set some parameters
#' growth <- function(t, x, params) { # A function returning a derivative
#'     with(params, {# r, K
#'         dx <- r*x*(1 - (x/K))
#'         return(list(dx))
#'    })
#'}
#'control <- list(times = 0:10, deltaT = 0.01)
#' simtimes <- seq(control$times[1], control$times[length(control$times)], by = control$deltaT) # for comparison between ode() and sde()
#' 
#' x.ode <- ode(.growth$xinit, simtimes, growth, .growth)
#' x.sde <- sde(.growth$xinit, control$times, growth, .growth, control)
#' 
#' plot(x.ode[, 1], x.ode[, 2], type = "l", xlab = "Time", ylab = "x")
#' lines(x.sde[, 1], x.sde[, 2], type = "l", col = 2)
#' 
#' 
#' .Brownian1D <- list(xinit = 0, mu = 0, sigma = 1e-1)
#' Brownian1D <- function(t, x, params) {
#'     with(params, {
#'         dx <- mu*x
#'         return(list(dx))
#'     })
#' }
#' 
#' .Brownian2D <- list(xinit = c(y = 0, z = 0), mu = 0, sigma = 1e-1, N = 1, nstatevars = 2)
#' Brownian2D <- function(t, x, params) {
#'     with(params, {
#'         y <- x[1:N]
#'         z <- x[(N+1):(2*N)]
#' 
#'         dy <- mu*y
#'         dz <- mu*z
#'         return(list(dy, dz))
#'     })
#' }
#' 
#' x <- sde(.Brownian1D$xinit, control$times, Brownian1D, .Brownian1D, control)
#' plot(x[, 1], x[, 2], type = "l", xlab = "Time", ylab = "x")
#' 
#' res <- sde(.Brownian2D$xinit, control$times, Brownian2D, .Brownian2D, control)
#' 
#' matplot(res[, 1], res[, 2:3], type = "l", lty = 1, col = 1, lwd = 1, xlab = "Time", ylab = "x")
#' 
#' xlim <- ylim <- c(-max(res[, 2:3]), max(res[, 2:3]))
#' plot(NULL, xlim = xlim, ylim = ylim, xlab = "x", ylab = "y")
#' abline(h = 0, lwd = 0.5)
#' abline(v = 0, lwd = 0.5)
#' lines(res[, 2], res[, 3])
#' points(res[1, 2], res[1, 3], col = 3, pch = 16, cex = 2) # starting point
#' points(res[nrow(res), 2], res[nrow(res), 3], col = 2, pch = 16, cex = 2) # stopping point
#'
#' ## Coupled SDEs on a network with an absorbing state.
#' 
#' library(igraph)
#'
#' g <- make_full_graph(4)
#' A <- as_adj(g, "both", sparse = FALSE)
#' N <- vcount(g)
#' model <- SIS
#' params <- c(.SIS, list(A = A))
#' control <- list(
#'   deltaT = 0.01, times = 0:50,
#'   absorbing.state = list(value = 0, which = "floor")
#' )
#' X <- sde(rep(0.01, N), control$times, model, params, control)
#' time_ev(X, ylim = c(-0.001, 0.001))
#' abline(h = 0, col = 2)
#' min(X) # 0
#' @export
sde <- function(initialvalue, times, func, parms = list(), control = list()) { # `parms` b/c deSolve
                                        # must have a noise strength
    stopifnot("sigma" %in% names(parms))
                                        # must have a Δt
    stopifnot("deltaT" %in% names(control))
                                        # must have adjacency matrix for coupled SDEs, but this doesn't enforce it
    if("A" %in% names(parms)) N <- nrow(parms$A) else N <- 1
    if("nstatevars" %in% names(parms)) nstatevars <- parms$nstatevars else nstatevars <- 1

    ntimesteps <- determine_ntimesteps(times, control$deltaT)
    W <- preallocate_noise(parms$sigma, ntimesteps, N*nstatevars)

    simtimes <- seq_len(ntimesteps)
    showtimes <- (simtimes - 1)*control$deltaT # Take a note of this: times arg is [minT, maxT). I think this makes sense: We start at 0, store our location at the initial point, then start moving. I want a trace since the beginning. Doesn't affect the output of solve_in_range(). 

    if(is.null(names(initialvalue))) varnames <- 1:N else varnames <- names(initialvalue)
    ##if(nstatevars > 1) {
        ##varnames <- 
    
    x <- initialvalue
    X <- matrix(0, nrow = ntimesteps, ncol = N*nstatevars)
    for(timestep in simtimes) {
        X[timestep, ] <- x
        x <- x +
                                        # This step doesn't use deSolve, only the deSolve format
                                        # i.e., functions which are written to be used with deSolve
                                        # the func() here is just a function that computes a next step
                                        # (though only tested with a derivative-computing function).
                                        # scaling for Δt for the derivative
            func(timestep, x, parms)[[1]]*control$deltaT + # from deSolve docs, derivative is [[1]]
                                        # scaling for Δt for the noise process
            W[timestep, ]*sqrt(control$deltaT)

        if("absorbing.state" %in% names(control)) {
            x <- enforce_absorbing_state(
                x,
                val = control$absorbing.state$value, which = control$absorbing.state$which)
        }
    }

    ##df <- as.data.frame(cbind(showtimes, X))
    res <- cbind(showtimes, X)
    ## colnames(df) <- c("time", varnames)
    colnames(res) <- c("time", varnames)
    ## return(df)
    return(res)
}

#' Handle an absorbing state in an SDE model
#'
#' Require that the output of sde() respect a value which should, in an ODE, be an absorbing state. Without this enforcement, the x_i may become larger or smaller than the absorbing state due to dynamical noise. Intended to be called from within sde() by using a named list within the control list. 
#' @param x The x_i returned by the deterministic part of the model plus the stochastic part
#' @param val The numeric value of the absorbing state (e.g., 0)
#' @param which Either "floor" or "ceiling
#' @returns A numeric vector, the x_i
#' @details Some models have absorbing states below or above which values should not stray. In deterministic simulations, the x_i stay where they should. In stochastic simulations, dynamical noise may result in x_i outside of the prescribed boundaries. This function enforces the absorbing state. The solution used is to reset the offending values by brute force. For example, if which = "floor", this is done: ifelse(x < val, val, x).
#'
#' It is not intended that this function is used on its own, but rather inside sde(). This is accomplished by means of a list with a particular name, "absorbing.state", and two named values, "value" and "which", e.g. control <- list(..., absorbing.state = list(value = 0, which = "floor")). See the examples.
#' @examples
#' library(graphics)
#' library(igraph)
#' library(localsolver)
#'
#' g <- make_full_graph(4)
#' A <- as_adj(g, "both", sparse = FALSE)
#' N <- vcount(g)
#' model <- SIS
#' params <- c(.SIS, list(A = A))
#' control <- list(
#'   deltaT = 0.01, times = 0:50,
#'   absorbing.state = list(value = 0, which = "floor")
#' )
#' X <- sde(rep(0.01, N), control$times, model, params, control)
#' time_ev(X, ylim = c(-0.001, 0.001))
#' abline(h = 0, col = 2)
#' min(X) # 0
#' 
#' X.wrong <- sde(rep(0.01, N), control$times, model, params, list(deltaT = 0.01))
#' min(X.wrong) # a small negative value
#'
#' dev.new(width = 14)
#' par(mfrow = c(1, 2))
#' time_ev(X, main = "Right", ylim = c(-1e-4, 1e-4))
#' abline(h = 0, col = 2)
#' time_ev(X.wrong, main = "Wrong", ylim = c(-1e-4, 1e-4))
#' abline(h = 0, col = 2)
#'
#' model <- mutualistic
#' params <- c(.mutualistic, list(A = A))
#' params$D <- 0.05
#' params$u <- -5
#' control$times <- 0:10
#' X <- sde(rep(10, N), control$times, model, params, control)
#' time_ev(X, ylim = c(-0.001, 0.001))
#'
#' model <- genereg
#' params <- c(.genereg, list(A = A))
#' params$D <- 0.25
#' control$times <- 0:50
#' X <- sde(rep(10, N), control$times, model, params, control)
#' time_ev(X, ylim = c(-1e-4, 1e-4))
enforce_absorbing_state <- function(x, val, which = c("floor", "ceiling")) {
    which <- match.arg(which)

    switch(
        which,
        floor = ifelse(x < val, val, x),
        ceiling = ifelse(x > val, val, x)
    )
}
