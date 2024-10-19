#' Functions to standardize simulations of ODE/SDEs
#'
#' Paired functions and lists of parameters to standardize using deSolve's ode() to compute ground truth simulations/solutions of several ODEs on various networks. 
#' @name ODEs
#' @param t An arbitrary time.
#' @param x The current state of the system. When N > 1, x is a vector with x_i representing the current state of the i'th node. 
#' @param params A list of model parameters.
#' @details Names without dots are functions which compute derivatives and are in deSolve's standard form. Each function can be used in a one-time way to compute the derivative of a system given an arbitrary time t, current state x, and parameters. More commonly, passed to deSolve's ode() or sdn's sde() as the model to be simulated.
#' Dot names are lists of standard model parameters. The adjacency matrix is also required for solutions on networks; concatenate it to the params list like `c(params, list(A = A))`. If analyzing a single variable, x should have length 1 and pass A = matrix(0).
#' @return A vector of derivatives
#' @examples
#' library(parallel)
#' ncores <- detectCores()-1
#' library(igraph)
#' library(deSolve)
#' library(sdn)
#' g <- sample_pa(50, m = 2, directed = FALSE, start.graph = make_full_graph(3))
#' N <- vcount(g)
#' A <- as_adj(g, "both")#, sparse = FALSE)
#' times <- 0:15
#' params <- c(.doublewell, list(A = A))
#' control <- list(times = times, ncores = ncores)
#' system.time(
#'   X <- solve_in_range(params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode")
#' ) # 63 seconds with 4 i3-5010U CPU @ 2.10 GHz
#' bifplot(X, params$Ds, TRUE, col = adjustcolor(1, 0.25))
NULL

#' @rdname ODEs
#' @export
.doublewell <- list(
    xinit.low = 1, xinit.high = 5,
    r = c(1, 3, 5),
    D = 0.05, Ds = seq(0, 1, length.out = 100),
    u = 0, us.up = seq(0, 5, length.out = 100), us.down = seq(0, -5, length.out = 100),
    sigma = 1e-2,
    basin.limit = 3 # separatrix
)

#' @rdname ODEs
#' @export
doublewell <- function(t, x, params) {
    with(params, {
        ## coupling <- D*rowSums(A*outer(rep(1, length(x)), x))
        coupling <- D*sapply(AL, function(vs) sum(x[vs]))
        dx <- -(x - r[1])*(x - r[2])*(x - r[3]) + coupling + u
        return(list(c(dx)))
    })
}

#' @rdname ODEs
#' @export
.SIS <- list(
    xinit.low = 0.001, xinit.high = 0.999,
    mu = 1,
    D = 0.05, Ds = seq(0, 1, length.out = 100),
    u = 0,
    sigma = 1e-3,#1e-5,
    basin.limit = 5e-3 # epidemic threshold, 5*sigma, clearly distinguishable from zero
)

#' @rdname ODEs
#' @export
SIS <- function(t, x, params) {
    with(params, {
        ## coupling <- D*rowSums(A*outer(1 - x, x))
        coupling <- D*mapply(function(xi, vs) sum((1 - xi)*x[vs]), x, AL)
        dx <- coupling - mu*x
        return(list(c(dx)))
    })
}

#' @rdname ODEs
#' @export
.genereg <- list( # Only the down direction makes sense for this model (going towards cell death)
    xinit.high = 2,
    B = 1, f = 1, h = 2,
    D = 1, Ds = seq(0, 1, length.out = 100),
    u = 0, us.down = seq(0, -1, length.out = 100),
    sigma = 1e-3,
    basin.limit = 5e-3 # cell death, 5*sigma, clearly distinguishable from zero
)

#' @rdname ODEs
#' @export
genereg <- function(t, x, params) {
    with(params, {
        ## coupling <- D*rowSums(A*outer(rep(1, length(x)), (x^h)/(1 + (x^h))))
        coupling <- D*sapply(AL, function(vs) sum(x[vs]^h/(1 + x[vs]^h)))
        dx <- ifelse( # handles situation where control parameter is u. If sde, use absorbing.state list in control
            x > 0,
            -B*(x^f) + coupling + u,
            0
        )
        return(list(c(dx)))
    })
}

#' @rdname ODEs
#' @export
.mutualistic <- list( # Up direction could be invasion, down collapse
    xinit.low = 0.001, xinit.high = 6, # to match Gao et al 2016
    B = 0.1, K = 5, C = 1, Dtilde = 5, E = 0.9, H = 0.1,
    D = 0.05, Ds = seq(0, 3, length.out = 100),
    u = 0, us.up = seq(0, 0.5, length.out = 100), us.down = seq(0, -5, length.out = 100),
    sigma = 1e-3,
    basin.limit = 1 # Allee constant (= C)
)

#' @rdname ODEs
#' @export
mutualistic <- function(t, x, params) {
    with(params, {
        ## coupling <- D*rowSums(A*(outer(x, x)/(Dtilde + outer(E*x, H*x, `+`))))
        coupling <- D*mapply(function(xi, vs) sum((xi*x[vs])/(Dtilde + E*xi + H*x[vs])), x, AL)
        dx <- ifelse(
            x > 0,
            B + x*(1 - (x/K))*((x/C) - 1) + coupling + u,
            0
        )
        return(list(c(dx)))
    })
}
