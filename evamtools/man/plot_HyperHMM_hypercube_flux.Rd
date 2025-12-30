\name{plot_HyperHMM_hypercube_flux}
\alias{plot_HyperHMM_hypercube_flux}
\title{Plot Hypercube Flux from HyperHMM}
\description{
    This function visualized the progression pathway inferred by HyperHMM as a DAG (directed acyclic graph).
    It represents the hypercube of states where nodes represent genetic states and edges represent transitions betweem them.
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

    \item{node.labels}{Logical. If \code{TRUE} (default), the binary or decimal names of the model will 
    be displayed.}

    \item{use.probability}{Logical.}
}