\name{plot_HyperHMM_sample_genots}
\alias{plot_HyperHMM_sample_genots}
\title{Plot of sampled genotypes}
\description{
Plot of sampled genotypes from cross-sectional data.
This function is adapted from internal plotting utilities and uses genotype
counts to display their counts. 
The colour of the plot can be customised.
}

\usage{
plot_HyperHMM_sample_genots(cpm_output)
}

\arguments{
\item{cpm_output}{Output from the cpm (result of running \code{evam()}.)}
}

\value{
    A base R plot showing sampled genotype counts.
}

\examples{
example_matrix <- matrix(c(0,0,1,0,1,1,0,0,0,1,0,1), byrow = TRUE, ncol = 3)
colnames(example_matrix) <- paste0("Gene", 1:3)
model_HyperHMM <- evam(example_matrix, methods = "HyperHMM")  
plot_HyperHMM_sample_genots(model_HyperHMM)
}

\references{
    \bold{evam} 

- Diaz-Uriarte, R & Herrera-Nieto, P. 2022. EvAM-Tools: tools for evolutionary accumulation and cancer 
progression models. Bioinformatics. url\{https://doi.org/10.1093/bioinformatics/btac710}
}