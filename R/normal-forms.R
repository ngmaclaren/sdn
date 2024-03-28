#' Saddle-node bifurcations
#'
#' @name saddlenode
#' @description The normal form of the saddle-node bifurcation.
#'
#' @param r The control parameter.
#'
#' @references ${1: https://en.wikipedia.org/wiki/Saddle-node_bifurcation}
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
#'
#' @param r The control parameter.
#'
#' @references ${1: https://en.wikipedia.org/wiki/Transcritical_bifurcation}
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
#' @description The normal form of the pitchfork bifurcation.
#' @param r The control parameter.
#' @references ${1: https://en.wikipedia.org/wiki/Pitchfork_bifurcation}
NULL
#' @rdname pitchfork
#' @export
.pitchfork <- list(
    xinit = 5, rs = seq(-5, 5, length.out = L), threshold = 0.01,
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
#' @details Hopf bifurcations are necessarily bivariate, requiring both an x and a y coordinate for a single node or isolated system. The implementation here follows from that used in deSolve: the 'x' variable is a vector with the x and y values concatenated. The Hopf() function itself splits the x and y values, computes the derivative, the concatenates and returns.
#' @references ${1: https://en.wikipedia.org/wiki/Hopf_bifurcation}
NULL
#' @rdname Hopf
#' @export
.Hopf <- list(
    xinit = c(y = 0.5, z = 0.5), nstatevars = 2,
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
