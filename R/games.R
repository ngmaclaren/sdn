## Follow Allen et al 2017 to simulate evolutionary games on arbitrary networks
## Death-birth updating
## Assumes undirected graphs, but I guess that could be changed?
##
## Required parameters:
#### η :: selection strength
#### b :: benefit in donation game
#### c :: cost in donation game
#### g :: network (not really a parameter: it's the space)
## Seed
#### x_i :: node states

## TODO
## Check with weights. No expectations, though? Just need to make sure the weights are being used.
## Roxygenize.

## assumes igraph
## should probably parallelize...

calc_rrate <- function(u, g, params) { # u is a igraph.vs, g is an igraph, eta is a numeric (scalar)
    if(!inherits(u, "igraph.vs")) u <- V(g)[u]

    nbs <- neighbors(g, u)
    neighbors_donating <- sum(V(g)[nbs]$state)

    w_ij <- E(g)[u %--% nbs]$weight
    s_i <- strength(g, u)
    p_ij <- w_ij/s_i

    with(params, 1 + eta*(-c*u$state + b*sum(p_ij*V(g)[nbs]$state)))
}
        
death_birth <- function(u, g, params) {
                                        # node was passed to this function, probably randomly selected but it
                                        # doesn't have to have been.
                                        # u is the dying node
    if(!inherits(u, "igraph.vs")) u <- V(g)[u]
    
    nbs <- neighbors(g, u)
    if(length(nbs) == 1) {
        v <- V(g)[nbs]
    } else {
        rrates <- sapply(nbs, calc_rrate, g, params)
        w_ij <- E(g)[u %--% nbs]$weight
        
        probs <- w_ij*rrates/sum(w_ij*rrates)
        v <- V(g)[sample(nbs, 1, prob = probs)]
    }

    return(v$state) # so, select the node, then call death_birth() on that node, possibly changing its state
}
    
## format: y, times, func, parms, method, ...
ede <- function(statevector, g, simtime, update, params, keep_history = FALSE) {
    stopifnot(length(statevector) == vcount(g))
    stopifnot("weight" %in% names(edge_attr(g)))

    V(g)$state <- statevector
    V(g)$strength <- strength(g)
    N <- vcount(g)
    
    if(keep_history) S <- matrix(0, nrow = simtime, ncol = N)
    

    for(st in seq(simtime)) {
        if(keep_history) S[st, ] <- V(g)$state
        if(all(V(g)$state == 0) | all(V(g)$state == 1)) break
        
        u <- V(g)[sample.int(N, 1)]
        V(g)[u]$state <- update(u, g, params) # death_birth
    }

    if(keep_history) {
        return(S)
    } else {
        return(V(g)$state)
    }
}
