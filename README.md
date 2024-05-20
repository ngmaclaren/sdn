# sdn

An R package for simulating dynamics on networks. The primary motivations for this package were: (1) simulating dynamical equations coupled through an adjacency matrix (or Laplacian matrix, but I haven't tried it yet), and (2) the need to switch back and forth between deterministic and stochastic simulations of the same dynamics. For ODEs, relies on the [deSolve](https://cran.r-project.org/package=deSolve) package. For SDEs, implements some of the ideas in the [sde](https://cran.r-project.org/package=sde) package and its accompanying [book](https://link.springer.com/book/10.1007/978-0-387-75839-8) with the goal of interoperability with deSolve. 

## Installation

For now, I recommend cloning this repository, then using the following shell commands to build and install the package (assuming your working directory is the one that contains the cloned sdn directory):

```sh
R CMD build sdn
R CMD INSTALL sdn_0.1-1.tar.gz
```

