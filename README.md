# localsolver

A simple package to help automate some differential equation solving tasks in R. For ODEs, relies on the [deSolve](https://cran.r-project.org/package=deSolve) package. For SDEs, implements some of the ideas in the [sde](https://cran.r-project.org/package=sde) package and its accompanying [book](https://link.springer.com/book/10.1007/978-0-387-75839-8) with the goal of interoperability with deSolve. 

## Installation

For now, I recommend cloning this repository, then using the following shell commands to build and install the package (assuming your working directory is the one that contains the cloned localsolver directory):

```sh
R CMD build localsolver
R CMD INSTALL localsolver_0.1-1.tar.gz
```

