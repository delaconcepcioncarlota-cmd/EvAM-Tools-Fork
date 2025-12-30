\name{plot_HyperHMM_hypercube_flux}
\alias{plot_HyperHMM_hypercube_flux}
\title{Adapted HyperTraPS plot for HyperHMM outputs}
\description{
    This function visualized the progression pathway inferred by HyperHMM as a DAG (directed acyclic graph).
    It represents the hypercube of states where nodes represent genetic states and edges represent transitions betweem them.
    This function has been adapted from the HyperHMM package, available at \url{https://github.com/StochasticBiology/hypercube-hmm.git}.
}

\usage{
    plot_HyperHMM_hypercube_flux(
        fitted.obj,
        thresh = 0.05,
        node.labels = TRUE
        use.probability = FALSE
    )
}

\arguments{
    \item{fitted.obj}{The output object directly returned by the \code{HyperHMM} function.
    It must contain the transition data frame (transitions) and the number of features (L).}

    \item{thresh}{A numeric velue (0.05 as default). Transitions with a Flux or Probability 
    below this threshold will be filtered out to simplify the visualization.}

    \item{node.labels}{Boolean. default TRUE. The binary or decimal names of the model will 
    be displayed or not.}

    \item{use.probability}{Boolean. default FALSE. It use probability instead of flux or not.}

    \value{
        A ggplot output
    }
}
\examples{
m = matrix(c(0,0,1,0,1,1), byrow=TRUE, ncol=3)
fitted = EvamTools(m, model = "HyperHMM")
plot_HyperHMM_hypercube_flux(fitted)
}

\references{
    \bold{HyperHMM} 

- Moen, M. T., & Johnston, I. G. (2022). Efficient inference of evolutionary and progressive dynamics 
on hypercubic transition graphs. bioRxiv, 2022.05.09.491130. \url https://doi.org/10.1101/2022.05.09.491130 

}
