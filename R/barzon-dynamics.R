#' Biochemical
#'
#' A simple interaction model.
#' @name biochemical
#' @param t An arbitrary time.
#' @param x The current state of the system. When N > 1, x is a vector with x_i representing the current state of the i'th node.
#' @param params A list of model parameters.
#' @details A simple model in deSolve's format.
#' @return A vector of derivatives.
#' @references Harush
#' @export
biochemical <- function(t, x, params) {
    with(params, {
        coupling <- R*mapply(function(xi, vs) sum(xi*x[vs]), x, AL)
        dx <- F - B*x - coupling
        return(list(c(dx)))
    })
}

#' @rdname biochemical
#' @export
.biochemical <- list(
    F = 1, B = 1, R = 1,
    xinit = 0.001
)  

#' Epidemic
#'
#' A susceptible-infected-susceptible (SIS) model on networks.
#' @name epidemic
#' @param t An arbitrary time.
#' @param x The current state of the system. When N > 1, x is a vector with x_i representing the current state of the i'th node.
#' @param params A list of model parameters.
#' @details A simple model in deSolve's format.
#' @return A vector of derivatives.
#' @references Harush
#' @export
epidemic <- function(t, x, params) {
    with(params, {
        coupling <- R*mapply(function(xi, vs) sum((1 - xi)*x[vs]), x, AL)
        dx <- -B*x + coupling
        return(list(c(dx)))
    })
}

#' @rdname epidemic
#' @export
.epidemic <- list(
    B = 1, R = 1,
    xinit = 0.001
)

#' Reduced mutualistic
#'
#' A reduced (in the sense of having few parameters) model of mutualistic dynamics.
#' @name mutual
#' @param t An arbitrary time.
#' @param x The current state of the system. When N > 1, x is a vector with x_i representing the current state of the i'th node.
#' @param params A list of model parameters.
#' @details A simple model in deSolve's format.
#' @return A vector of derivatives.
#' @references Harush
#' @export
mutual <- function(t, x, params) {
    with(params, {
        coupling <- R*mapply(function(xi, vs) sum(xi*(x[vs]^b/(1 + x[vs]^b))), x, AL)
        dx <- B*x*(1 - x) + coupling
        return(list(c(dx)))
    })
}

#' @rdname mutual
#' @export
.mutual <- list(
    B = 1, b = 2, R = 1,
    xinit = 0.001
)
