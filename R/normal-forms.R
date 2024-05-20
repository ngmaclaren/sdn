#' Saddle-node bifurcations
#'
#' @name saddlenode
#' @description The normal form of the saddle-node bifurcation.
#' @param r The control parameter.
#' @references https://en.wikipedia.org/wiki/Saddle-node_bifurcation
#' @examples
#' library(deSolve)
#' library(localsolver)
#' rs <- .saddlenode$rs
#' x <- solve_in_range(rs, "r", saddlenode, .saddlenode$xinit, .saddlenode)
#' xlim <- range(rs) + c(-0.1, 0.1)
#' ylim <- c(min(x), max(c(max(x), 0.1)))
#' plot(rs, x, type = "l", xlab = "r", ylab = "x", col = 2, lwd = 2,
#' xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")
#' abline(h = 0)
#' abline(v = 0)
NULL

#' @rdname saddlenode
#' @export
.saddlenode <- list(
    xinit = -5, rs = seq(-2.5, 0, length.out = 100),
    sigma = 1e-2
)
#' @rdname saddlenode
#' @export
saddlenode <- function(t, x, params) {
    with(params, {
        dx <- r + x^2
        return(list(c(dx)))
    })
}

#' Transcritical bifurcations
#'
#' @name transcritical
#' @description The normal form of the transcritical bifurcation
#' @param r The control parameter.
#' @references https://en.wikipedia.org/wiki/Transcritical_bifurcation
#' @examples
#' library(deSolve)
#' library(localsolver)
#' rs <- .transcritical$rs
#' x <- solve_in_range(rs, "r", transcritical, .transcritical$xinit, .transcritical)
#' plot(rs, x, type = "l", xlab = "r", ylab = "x", col = 2, lwd = 2, xaxs = "i")
#' abline(h = 0)
#' abline(v = 0)
NULL

#' @rdname transcritical
#' @export
.transcritical <- list(
    xinit = 5, rs = seq(-5, 5, length.out = 100), threshold = 0.01,
    sigma = 1e-2
)

#' @rdname transcritical
#' @export
transcritical <- function(t, x, params) {
    with(params, {
        dx <- r*x - x^2
        return(list(c(dx)))
    })
}

#' Pitchfork bifurcations
#' @name pitchfork
#' @description The normal form of the pitchfork bifurcation. This function is the supercritical type.
#' @param r The control parameter.
#' @references https://en.wikipedia.org/wiki/Pitchfork_bifurcation
#' @examples
#' library(graphics)
#' library(deSolve)
#' library(localsolver)
#' rs <- .pitchfork$rs
#' x_upper <- solve_in_range(rs, "r", pitchfork, .pitchfork$xinit, .pitchfork)
#' x_lower <- solve_in_range(rs, "r", pitchfork, -.pitchfork$xinit, .pitchfork)
#' ylim <- range(c(x_upper, x_lower))
#' plot(rs, x_upper, type = "l", xlab = "r", ylab = "x", col = 2, lwd = 2, xaxs = "i", ylim = ylim)
#' lines(rs, x_lower, col = 3, lwd = 2)
#' abline(h = 0)
#' abline(v = 0)
#'
#' ## Note that increasing the simulation time will increase the resolution on the bifurcation point
#' control <- list(times = 0:100)
#' x_upper <- solve_in_range(rs, "r", pitchfork, .pitchfork$xinit, .pitchfork, control)
#' x_lower <- solve_in_range(rs, "r", pitchfork, -.pitchfork$xinit, .pitchfork, control)
#' ylim <- range(c(x_upper, x_lower))
#' plot(rs, x_upper, type = "l", xlab = "r", ylab = "x", col = 2, lwd = 2, xaxs = "i", ylim = ylim)
#' lines(rs, x_lower, col = 3, lwd = 2)
#' abline(h = 0)
#' abline(v = 0)
NULL
#' @rdname pitchfork
#' @export
.pitchfork <- list(
    xinit = 5, rs = seq(-5, 5, length.out = 100), threshold = 0.01,
    sigma = 1e-2
)
#' @rdname pitchfork
#' @export
pitchfork <- function(t, x, params) {
    with(params, {
        dx <-r*x - x^3
        return(list(c(dx)))

    })
}

#' Hopf bifurcations
#' @name Hopf
#' @description The normal form of the Hopf bifurcation in Euclidean coordinates and without complex numbers.
#' @param r One of two possible control parameters, a.k.a. mu
#' @param w One of two possible control parameters, a.k.a. omega
#' @details Hopf bifurcations are necessarily bivariate, requiring both an x and a y coordinate for a single node or isolated system. The implementation here follows from that used in deSolve: the 'x' variable is a vector with the x and y values concatenated. The Hopf() function itself splits the x and y values, computes the derivative, then concatenates and returns.
#' @references https://en.wikipedia.org/wiki/Hopf_bifurcation
#' @examples
#' library(deSolve)
#' library(localsolver)
#' rs <- c(-1, -0.1, -0.01, -0.001, 0.001, 0.01, 0.1, 1)
#' control <- list(times = seq(0, 100, by = 0.1))
#' result <- lapply(rs, function(r) {
#'   ode(.Hopf$xinit, control$times, Hopf, c(list(r = r), .Hopf))
#' })
#' plotit <- function(res) {
#'   start <- res[1, c("y", "z")]
#'   stop <- res[nrow(res), c("y", "z")]
#'   xlim = c(-1, 1)
#'   ylim = c(-1, 1)
#'   plot(
#'     res[, c("y", "z")], type = "l", col = 8, lwd = 1,
#'     xlab = "x", ylab = "y", xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i"
#'   )
#'   points(start["y"], start["z"], pch = 16, col = 3, cex = 3)
#'   points(stop["y"], stop["z"], pch = 16, col = 2, cex = 3)
#' }
#' dev.new(height = 8, width = 16)
#' par(mfrow = c(2, 4))
#' for(i in seq_along(result)) {
#'   plotit(result[[i]])
#'   title(main = paste("r =", rs[i]))
#' }
NULL
#' @rdname Hopf
#' @export
.Hopf <- list(
    xinit = c(y = 0.5, z = 0.5), nstatevars = 2, N = 1,
    w = 1,
    rs = seq(-1, 1, length.out = 100),
    threshold = list(zero = 0.01, limitcycle = 0.9),
    sigma = 1e-2
)
#' @rdname Hopf
#' @export
Hopf <- function(t, x, params) {
    ## r is μ, w is ω
    with(params, {
        y <- x[1:N] # if only one node, N will be 1 and 1:1 is 1
        z <- x[(N+1):(2*N)] # if only one node, N+1 = 2 and 2*N = 2, z is a vector with length 2 ###?

        dy <- (r - y^2 - z^2)*y - w*z # ω = 1
        dz <- (r - y^2 - z^2)*z + w*y # ω = 1
        return(list(c(dy, dz)))
    })
}
