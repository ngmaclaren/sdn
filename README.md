# sdn

An R package for simulating dynamics on networks. The primary motivations for this package were: (1) simulating dynamical equations coupled through an adjacency matrix (or Laplacian matrix, but I haven't tried it yet), and (2) the need to switch back and forth between deterministic and stochastic simulations of the same dynamics. For ODEs, relies on the [deSolve](https://cran.r-project.org/package=deSolve) package. For SDEs, implements some of the ideas in the [sde](https://cran.r-project.org/package=sde) package and its accompanying [book](https://link.springer.com/book/10.1007/978-0-387-75839-8) with the goal of interoperability with deSolve. 

Version 0.2-1 introduces a breaking change: the default models, which I used in several projects, have been switched from adjacency matrix to adjacency list computation. The new version should be much faster. 

## Models

Recently updated to include models from [Barzon et al. 2024](https://github.com/gbarzon/jacobian_geometry), but check parameter values and ranges to confirm desired behavior. Code for new models is in `./R/barzon-dynamics.R`. Other models can be found in `./R/models.R` and `./R/normal-forms.R`. 

## Installation

For now, I recommend cloning this repository, then using the following shell commands to build and install the package (assuming your working directory is the one that contains the cloned sdn directory):

```sh
R CMD build sdn
R CMD INSTALL sdn_0.2-1.tar.gz
```

