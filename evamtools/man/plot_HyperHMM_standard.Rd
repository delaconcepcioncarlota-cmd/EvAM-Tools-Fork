\name{plot_HyperHMM_standard}
\alias{plot_HyperHMM_standard}
\title{Panel plot for model output}
\description{
    This function has been adapted from the HyperHMM package, available at \url{https://github.com/StochasticBiology/hypercube-hmm.git}.
    It represents a panel with three graphs. The hypercube_flux plot, the bubble plot and the PFG.
}

\usage{
plot_HyperHMM_standard(cpm_output, legacy = FALSE, label = "")
}

\arguments{
\item{cpm_output}{Output from the cpm. Only the \code{stats} component of \code{HyperHMM_primary_output} is used.}

\item{legacy}{redundant, unused}

\item{label}{(character) optional title for plot}
}
\value{
A ggplot output
}

\examples{
m = matrix(c(0,0,1,0,1,1), byrow=TRUE, ncol=3)
out = evam(m, methods = "HyperHMM")
plot_HyperHMM_standard(out)
}

\references{
    \bold{HyperHMM} 

- Moen, M. T., & Johnston, I. G. (2022). Efficient inference of evolutionary and progressive dynamics 
on hypercubic transition graphs. bioRxiv, 2022.05.09.491130. \url{https://doi.org/10.1101/2022.05.09.491130}

}