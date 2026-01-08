\name{plot_HyperHMM_as_evam_trans_mat}
\alias{plot_HyperHMM_as_evam_trans_mat}
\title{Plot genotype transition matrix from a HyperHMM model}
\description{Visualize the evolutionary pathways estimated by aHyperHMM model}

\usage{plot_HyperHMM_as_evam_trans_mat(
  cpm_output,
  label_type = "genotype",
  fixed_vertex_size = FALSE,
  top_paths = NULL
)
}

\arguments{
    \item{cpm_output}{Output from the cpm}
    \item{label_type}{Type of label to show. One of:
        \describe{
        \item{genotype: }{Displays all genes mutated}
        \item{acquisition: }{Only displays the last gente mutated}
    \item{fixed_vertex_size}{Boolean. If TRUE, all nodes with have the same
        size; otherwise, scale them proportional to frequencies of observed data.}
    \item{top_paths}{Number of most relevant paths to plot. Default NULL
        will plot all paths. See below, Description, for details about
        relevance. With many genes, and particularly for MHN (or other methods,
        when there are no restrictions in the order of accumulation of
        mutations), using NULL can lead to a very long time to plot.}

  }
}
}

\value{
    A igraph object
}

\examples{
    dB_c1 <- matrix(
    c(
        rep(c(1, 0, 0, 0, 0), 30) #A
    , rep(c(0, 0, 1, 0, 0), 30) #C
    , rep(c(1, 1, 0, 0, 0), 20) #AB
    , rep(c(0, 0, 1, 1, 0), 20) #CD
    , rep(c(1, 1, 1, 0, 0), 10) #ABC
    , rep(c(1, 0, 1, 1, 0), 10) #ACD
    , rep(c(1, 1, 0, 0, 1), 10) #ABE
    , rep(c(0, 0, 1, 1, 1), 10) #CDE
    , rep(c(1, 1, 1, 0, 1), 10) #ABCE
    , rep(c(1, 0, 1, 1, 1), 10) #ACDE
    , rep(c(1, 1, 1, 1, 0), 5) # ABCD
    , rep(c(0, 0, 0, 0, 0), 1) # WT
    ), ncol = 5, byrow = TRUE
    )
    colnames(dB_c1) <- LETTERS[1:5]
    out_h_2 <- evam(dB_c1,
                methods = c("HyperHMM"),
                hyper_hmm_opts = list(length = 2,
                                        model = 2,
                                        seed = -1))
    plot_HyperHMM_as_evam_trans_mat(out_h_2, top_paths = 4)
}

\references{
    \bold{evam} 

- Diaz-Uriarte, R & Herrera-Nieto, P. 2022. EvAM-Tools: tools for evolutionary accumulation and cancer 
progression models. Bioinformatics. url\{https://doi.org/10.1093/bioinformatics/btac710}
}