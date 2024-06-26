% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normal-forms.R
\docType{data}
\name{pitchfork}
\alias{pitchfork}
\alias{.pitchfork}
\title{Pitchfork bifurcations}
\format{
An object of class \code{list} of length 4.
}
\usage{
.pitchfork

pitchfork(t, x, params)
}
\arguments{
\item{r}{The control parameter.}
}
\description{
The normal form of the pitchfork bifurcation. This function is the supercritical type.
}
\examples{
library(graphics)
library(deSolve)
library(sdn)
rs <- .pitchfork$rs
x_upper <- solve_in_range(rs, "r", pitchfork, .pitchfork$xinit, .pitchfork)
x_lower <- solve_in_range(rs, "r", pitchfork, -.pitchfork$xinit, .pitchfork)
ylim <- range(c(x_upper, x_lower))
plot(rs, x_upper, type = "l", xlab = "r", ylab = "x", col = 2, lwd = 2, xaxs = "i", ylim = ylim)
lines(rs, x_lower, col = 3, lwd = 2)
abline(h = 0)
abline(v = 0)

## Note that increasing the simulation time will increase the resolution on the bifurcation point
control <- list(times = 0:100)
x_upper <- solve_in_range(rs, "r", pitchfork, .pitchfork$xinit, .pitchfork, control)
x_lower <- solve_in_range(rs, "r", pitchfork, -.pitchfork$xinit, .pitchfork, control)
ylim <- range(c(x_upper, x_lower))
plot(rs, x_upper, type = "l", xlab = "r", ylab = "x", col = 2, lwd = 2, xaxs = "i", ylim = ylim)
lines(rs, x_lower, col = 3, lwd = 2)
abline(h = 0)
abline(v = 0)
}
\references{
https://en.wikipedia.org/wiki/Pitchfork_bifurcation
}
\keyword{datasets}
