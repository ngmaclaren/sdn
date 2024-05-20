#' Plot the time evolution of a system
#' @description Given the output of deSolve::ode() or localsolver::sde(), plot the time evolution of the system.
#' @param X The output of ode() or sde()
#' @param showmean Should the mean state of potentially many variables in X be drawn?
#' @param ... Passed to matplot() and, if showmean = TRUE, lines()
#' @details This is meant as a quick diagnostic plots. At least one consequence of this is that `col` is hard-coded inside the function...
#' @examples
#' library(graphics)
#' library(deSolve)
#' library(localsolver)
#'
#' control <- list(times = seq(0, 20, by = 0.1))
#' x <- ode(.saddlenode$xinit, control$times, saddlenode, list(r = 0))
#' time_ev(x)
#' abline(h = 0, lty = 2)
#'
#' library(igraph)
#' g <- sample_pa(25, m = 2, directed = FALSE, start.graph = make_full_graph(3))
#' A <- as_adj(g, "both", sparse = FALSE)
#' N <- vcount(g)
#' params <- c(list(A = A), .doublewell)
#' control <- list(times = seq(0, 1, length.out = 100))
#' X <- ode(rep(.doublewell$xinit.low, N), control$times, doublewell, params)
#' time_ev(X, showmean = TRUE)
#' @export
time_ev <- function(X, showmean = FALSE, ...) {
    with(
        list(...), {
            matplot(X[, 1], X[, -1], type = "l", lty = 1, xlab = "t", ylab = "x", col = 1, ...)
            if(showmean) lines(X[, 1], rowMeans(X[, -1]), col = 2, lwd = 2, ...)
        }
    )
}

#' Plot the bifurcation diagram
#' @description Given the output of localsolver::solve_in_range(), plot the bifurcation diagram of the system as a function of a control parameter.
#' @param X The output of solve_in_range(), a matrix
#' @param cparam The control parameter. It should be a sequence.
#' @param showmean Should the mean state of potentially many variables in X be drawn?
#' @param ... Passed to matplot(); col = 2 is used for the mean state if showmean = TRUE
#' @details This is meant as a quick diagnostic plot.
#' @examples
#' library(graphics)
#' library(deSolve)
#' library(localsolver)
#' rs <- .saddlenode$rs
#' x <- solve_in_range(rs, "r", saddlenode, .saddlenode$xinit, .saddlenode)
#' bifplot(x, rs)
#'
#' library(igraph)
#' g <- sample_pa(25, m = 2, directed = FALSE, start.graph = make_full_graph(3))
#' A <- as_adj(g, "both", sparse = FALSE)
#' N <- vcount(g)
#' params <- c(list(A = A), .doublewell)
#' control <- list(times = seq(0, 1, length.out = 100))
#' X <- solve_in_range(params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control)
#' bifplot(X, params$Ds, showmean = TRUE, col = 1)
#' @export
bifplot <- function(X, cparam, showmean = FALSE, ...) {
                                        # Here X is the output of solve_in_range()
                                        # cparam should be a sequence
    with(list(...), {
        matplot(cparam, X, type = "l", lty = 1, xlab = "Control parameter", ylab = "x", ...)
    })
    if(showmean) lines(cparam, rowMeans(X), col = 2, lwd = 2)
}
