\name{plot_HyperHMM_pfg}
\alias{plot_HyperHMM_pfg}
\title{Plot graph of ordered pair acquisitions}
\description{
Plot graph of ordered pair acquisitions. Visualizes the Pathward Feature Graph (PFG), which represents
the preferred order of feature (mutation) acquisition. 
This function has been adapted from the HyperHMM package, available at \url{https://github.com/StochasticBiology/hypercube-hmm.git}.
}

\usage{
plot_HyperHMM_pfg(fitted.obj, pfg.layout = "matrix", curvature = 1)
}
\arguments{
\item{fitted.obj}{The output object directly returned by the \code{HyperHMM} function.
    It must contain the transition data frame (transitions) and the number of features (L).}

\item{pfg.layout}{Character string specifying the graph layout algorithm:
\itemize{
  \item \code{"matrix"}: (Default) Arranges nodes ina grid-like structure.
  \item \code{"tree"}: Attemps to arrange the features in a hierarchical tree-like structure.
  \item \code{"other"}: Uses a standard force-directed layout that prioritizes equal edge length.
}}

\item{curvature}{(double, default 1) curvature of edges}
}
\value{
A ggplot output
}

\examples{
m = matrix(c(0,0,1,0,1,1), byrow=TRUE, ncol=3)
fitted = evam(m, methods = "HyperHMM")
plot_HyperHMM_pfg(fitted$HyperHMM_primary_output)
}

\references{
    \bold{HyperHMM} 

- Moen, M. T., & Johnston, I. G. (2022). Efficient inference of evolutionary and progressive dynamics 
on hypercubic transition graphs. bioRxiv, 2022.05.09.491130. \url{https://doi.org/10.1101/2022.05.09.491130}

}
