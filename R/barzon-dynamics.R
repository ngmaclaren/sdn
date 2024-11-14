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

#' Population
#'
#' A simple population model.
#' @name population
#' @param t An arbitrary time.
#' @param x The current state of the system. When N > 1, x is a vector with x_i representing the current state of the i'th node.
#' @param params A list of model parameters.
#' @details A simple model in deSolve's format.
#' @return A vector of derivatives.
#' @references Harush
#' @export
population <- function(t, x, params) {
    with(params, {
        coupling <- R*sapply(AL, function(vs) sum(x[vs]^a))
        dx <- -B*x^b + coupling
        return(list(c(dx)))
    })
}

#' @rdname population
#' @export
.population <- list(
    R = 1, B = 1, b = 3, a = 2,
    xinit = 0.001
)

#' Neuronal
#'
#' A simple model of neuron activity(?).
#' @name neuronal
#' @param t An arbitrary time.
#' @param x The current state of the system. When N > 1, x is a vector with x_i representing the current state of the i'th node.
#' @param params A list of model parameters.
#' @details A simple model in deSolve's format.
#' @return A vector of derivatives.
#' @references Hens
#' @export
neuronal <- function(t, x, params) {
    with(params, {
        coupling <- R*sapply(AL, function(vs) sum(tanh(x[vs])))
        dx <- -B*x + C*tanh(x) + coupling
        return(list(c(dx)))
    })
}

#' @rdname neuronal
#' @export
.neuronal <- list(
    B = 2, C = 2.5, R = 1,
    xinit = 0.001
)

#' Noisy voter
#'
#' A simple noisy voter model.
#' @name voter
#' @param t An arbitrary time.
#' @param x The current state of the system. When N > 1, x is a vector with x_i representing the current state of the i'th node.
#' @param params A list of model parameters.
#' @details A simple model in deSolve's format. Use Atilde in .voter because A is the name of the adjacency matrix in sdn(). 
#' @return A vector of derivatives.
#' @references Carro
#' @export
voter <- function(t, x, params) {
    with(params, {
        coupling <- (C/k)*sapply(AL, function(vs) sum(x[vs]))
        dx <- Atilde - B*x + coupling
        return(list(c(dx)))
    })
}

#' @rdname voter
#' @export
.voter <- list(
    Atilde = 1, B = 1, C = 1,
    k = "supply the degree vector",
    xinit = 0.001
)

#' Regulatory
#'
#' A simple model of a regulatory system.
#' @param t An arbitrary time.
#' @param x The current state of the system. When N > 1, x is a vector with x_i representing the current state of the i'th node.
#' @param params A list of model parameters.
#' @details A simple model in deSolve's format.
#' @return A vector of derivatives.
#' @references Harush
#' @export
regulatory <- function(t, x, params) {
    with(params, {
        coupling <- R*sapply(AL, function(vs) sum(x[vs]^h/(1 + x[vs]^h)))
        dx <- -B*x^a + coupling
        return(list(c(dx)))
    })
}

#' @rdname regulatory
#' @export
.regulatory <- list(
    a = 1, h = 2, R = 1, B = 1,
    xinit = 2
)

#' Synchronization
#'
#' A simple Kuramoto model.
#' @param t An arbitrary time.
#' @param x The current state of the system. When N > 1, x is a vector with x_i representing the current state of the i'th node.
#' @param params A list of model parameters.
#' @details A simple model in deSolve's format. Make sure to change the default values of omega and xinit in .synchronization to observe non-constant behavior.
#' @return A vector of derivatives.
#' @references Barzon only?
#' @export
synchronization <- function(t, x, params) {
    with(params, {
        coupling <- R*mapply(function(xi, vs) sum(sin(x[vs] - xi)), x, AL)
        dx <- omega + coupling
        return(list(c(dx)))
    })
}

#' @rdname synchronization
#' @export
.synchronization <- list(
    omega = 1, R = 1,
    xinit = 1
)
