\name{plot_HyperHMM_bubbles}
\alias{plot_HyperHMM_bubbles}
\title{Bubble plot from fitted model}
\description{
Bubble plot from fitted model
This function has been adapted from the HyperHMM package, available at \url{https://github.com/StochasticBiology/hypercube-hmm.git}.
}

\usage{
plot_HyperHMM_bubbles(fitted.obj, labels = NULL, formatted = FALSE)
}
\arguments{
\item{fitted.obj}{The output object directly returned by the \code{HyperHMM} function.
    It must contain the transition data frame (transitions) and the number of features (L).}

\item{labels}{(vector of characters) optional set of labels for features in plot}

\item{formatted}{(Boolean, default FALSE) redundant, unused}
}
\value{
    A ggplot output
}

\examples{
m = matrix(c(0,0,1,0,1,1), byrow=TRUE, ncol=3)
fitted = evam(m, model = "HyperHMM")
plot_HyperHMM_bubbles(fitted$HyperHMM_primary_output)
}

\references{
    \bold{HyperHMM} 

- Moen, M. T., & Johnston, I. G. (2022). Efficient inference of evolutionary and progressive dynamics
on hypercubic transition graphs. bioRxiv, 2022.05.09.491130. \url{https://doi.org/10.1101/2022.05.09.491130}

}

