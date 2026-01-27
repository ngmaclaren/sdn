## Follow Allen et al 2017 to simulate evolutionary games on arbitrary networks
## Death-birth updating
## Assumes undirected graphs, but I guess that could be changed? (Allen et al assume undirected.)
##
## Required parameters:
#### η :: selection strength
#### b :: benefit in donation game
#### c :: cost in donation game
#### g :: network (not really a parameter: it's the space)
## Seed
#### x_i :: node states

## TODO
## [ ] Check with weights. No expectations, though? Just need to make sure the weights are being used.
## [ ] Roxygenize. Start a template with C-c C-o C-o

## assumes igraph
## should probably parallelize...

##' Calculate the reproduction rate
##'
##' @description Calculate the reproduction rate of a node. 
##' @param u The node. Can be an integer node index or an igraph.vs.
##' @param g The network, an igraph object.
##' @param params A list of parameters. Must include eta, b, and c.
##' @return A numeric vector of length 1.
##' @details The `params` list should include, at a minimum: eta, the selection strength; b, the benefit in the donation game; and c, the cost in the donation game. The essential computation is `1 + eta*(-c*u$state + b*sum(p_ij*V(g)[nbs]$state))`.
##' @references Allen et al. 2017, doi:10.1038/nature21723.
##' @export
calc_reprate <- function(u, g, params) { # u is a igraph.vs, g is an igraph, eta is a numeric (scalar)
    if(!inherits(u, "igraph.vs")) u <- V(g)[u]

    nbs <- neighbors(g, u)
    neighbors_donating <- sum(V(g)[nbs]$state)

    w_ij <- E(g)[u %--% nbs]$weight
    s_i <- strength(g, u)
    p_ij <- w_ij/s_i

    with(params, 1 + eta*(-c*u$state + b*sum(p_ij*V(g)[nbs]$state)))
}
        
##' Death-birth updating
##'
##' Performs the death-birth updating rule on a network.
##' @param u The node. Can be an integer node index or an igraph.vs.
##' @param g The network, an igraph object.
##' @param params A list of parameters. Must include eta, b, and c.
##' @return An integer vector of length 1, a node state.
##' @references Allen et al. 2017, doi:10.1038/nature21723.
##' @details The `params` list should include, at a minimum: eta, the selection strength; b, the benefit in the donation game; and c, the cost in the donation game.
##'
##' This function takes, as its main input, a selected node, u. Ordinarily this will be a randomly selected node, but the function doesn't require it. Then, a neighbor v of node u is selected at random, proportional to the reproductive rates of u's neighbors. This function returns the state of v, which should ordinarily be assigned to u, completing the death-birth cycle.
##' @export
death_birth <- function(u, g, params) {
                                        # node was passed to this function, probably randomly selected but it
                                        # doesn't have to have been.
                                        # u is the dying node
    if(!inherits(u, "igraph.vs")) u <- V(g)[u]
    
    nbs <- neighbors(g, u)
    if(length(nbs) == 1) {
        v <- V(g)[nbs]
    } else {
        ## this right here is actually the donation game (including calc_reprate())
        rrates <- sapply(nbs, calc_reprate, g, params)
        w_ij <- E(g)[u %--% nbs]$weight
        
        probs <- w_ij*rrates/sum(w_ij*rrates)
        v <- V(g)[sample(nbs, 1, prob = probs)]
    }

    return(v$state) # so, select the node, then call death_birth() on that node, possibly changing its state
}
    
## format: y, times, func, parms, method, ...

##' Evolutionary games on networks
##'
##' The goal of this function is to perform arbitrary games on a network of N nodes. 
##' @param statevector A vector of length N. Should be all 0 or 1.
##' @param g A network, should be an igraph object.
##' @param simtime Number of simulation steps, an integer.
##' @param update The update rule, a function.
##' @param params A list of parameters. Must include eta, b, and c.
##' @param keep_history Logical.
##' @return The final state vector of the network, or, if keep_history=TRUE, a matrix with `simtime` rows and `length(statevector)` columns.
##' @references Allen et al. 2017, doi:10.1038/nature21723.
##' @details The `params` list should include, at a minimum: eta, the selection strength; b, the benefit in the donation game; and c, the cost in the donation game.
##'
##' Given a vector of node states `statevector`, a network `g`, the number of simulation iterations to perform `simtime`, an update rule `update`, and a list of parameters `params`, simulates an evolutionary game following the methods described in Allen et al (2017). The simulation will stop when `simtime` is reached or if the state vector is either all 0 or all 1. 
##'
##' The name for this function follows the pattern ode (ordinary differential equations), sde (stochastic differential equations), ede (evolutionary...), but the 'de' doesn't mean anything. 
##' @export
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
