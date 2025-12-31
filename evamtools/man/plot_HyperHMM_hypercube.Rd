\name{plot_HyperHMM_hypercube}
\alias{plot_HyperHMM_hypercube}
\title{Hypercube plot from fitted model}
\description{
This function has been adapted from the HyperHMM package, available at \url{https://github.com/StochasticBiology/hypercube-hmm.git}.
It represents each possible path.
redundant? older less flexible version -- will be sparsely commented
}
\usage{
plot_HyperHMM_hypercube(
  fitted.obj,
  use.width = TRUE,
  duplicate.offset = 0,
  lab.size = 3,
  p.size = 1,
  node.labels = TRUE,
  seg.labels = TRUE,
  threshold = 0,
  break.redundancy = 0,
  rotate.phi = FALSE
)
}
\arguments{
\item{fitted.obj}{The output object directly returned by the \code{HyperHMM} function.}

\item{use.width}{(Boolean, default TRUE) use line width to display edge weights?}

\item{duplicate.offset}{(double, default 0) vertical offset for nodes in identical positions}

\item{lab.size}{(double, default 3) size for edge labels}

\item{p.size}{(double, default 1) point size}

\item{node.labels}{(Boolean, default TRUE) node labels, yes or no?}

\item{seg.labels}{(Boolean, default TRUE) line segment labels, yes or no?}

\item{threshold}{(double, default 0) ignore edges under a threshold in the hypercube plot}

\item{break.redundancy}{redundant, unused}

\item{rotate.phi}{(Boolean, default FALSE) rotate states out of the page (in case of trajectories bunched up near the top/bottom)}
}
\value{
  A ggplot output
}

\examples{
m = matrix(c(0,0,1,0,1,1), byrow=TRUE, ncol=3)
fitted = evam(m, methods = "HyperHMM")
plot_HyperHMM_hypercube(fitted$HyperHMM_primary_output)
}

\references{
    \bold{HyperHMM} 

- Moen, M. T., & Johnston, I. G. (2022). Efficient inference of evolutionary and progressive dynamics 
on hypercubic transition graphs. bioRxiv, 2022.05.09.491130. \url{https://doi.org/10.1101/2022.05.09.491130}

}