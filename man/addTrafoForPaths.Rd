% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toolsSeverin.R
\name{addTrafoForPaths}
\alias{addTrafoForPaths}
\title{transform parametersets from the profiles for path plotting}
\usage{
addTrafoForPaths(profs, trafo, rescale = c("lin", "log", "log10", "log2")[1])
}
\arguments{
\item{profs}{parframe with the profiles, as returned from the \link{dMod::profile} function}

\item{trafo}{parameter transformation for the steady states, as returned by \code{P(steadystateTrafo)}. Currently no ther formulation is supported.}

\item{rescale}{character, default \code{"lin"} (no rescaling). The rescaling of the transformed parameters to the model scale, can be \code{"lin"}, \code{"log"}, \code{"log10"} or \code{"log2"}.}
}
\value{
\code{parframe} of the input \code{profs} with the added columns of \code{trafo} applied to the parameters.
}
\description{
while using non-trivial steady states, parameters can couple trough the steady states. This function applies the transformation to the parameter sets from the profiles, to account for that.
}
