% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solvers.R
\name{determine_ntimesteps}
\alias{determine_ntimesteps}
\title{Determine the number of time steps for a simulation}
\usage{
determine_ntimesteps(times = numeric(), maxT = 1, minT = 0, deltaT = 0.01)
}
\arguments{
\item{times}{a sequence of time steps, in user units}

\item{maxT}{The last time step (in user units)}

\item{minT}{The first time step (in user units)}

\item{deltaT}{The size of the timestep for simulation}
}
\value{
A scalar
}
\description{
Uses the desired time span and the simulation time step to compute how many time steps are required for the simulation.
}
\examples{
determine_ntimesteps(times = 0:10, deltaT = 0.01)
}
