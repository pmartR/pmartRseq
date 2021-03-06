% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_env.R
\name{mod_env}
\alias{mod_env}
\title{Function to calculate module eigen vectors and correlate against environmental variables}
\usage{
mod_env(omicsData, modData, envVars, pca.method = "svd",
  cor.method = "spearman", use = "pairwise", padjust = "BH")
}
\arguments{
\item{omicsData}{An object of the class 'seqData' usually created by \code{\link{as.seqData}}. omicsData$f_data must contain environmental variables to compare against.}

\item{modData}{An object of class 'modData', created by \code{\link{detect_modules}}, if want to colour by modules.}

\item{envVars}{A character vector with the names of the environmental variables to compare against. Must all be column names in omicsData$f_data. If NULL, will not run correlation with environmental variables, but will still calculate pca for modules.}

\item{pca.method}{A string specifying the pca method to use. Default is 'svd'. Other options include 'nipals', 'rnipals', 'bpca', 'ppca', 'svdImpute', 'robusePca', 'nlpca', 'llsImpute', and 'llsImputeAll'.}

\item{cor.method}{A string specifying the correlation method to use. Default is 'spearman', which is nonparametric and doesn't make any assumptions about the data. Other options include 'kendall' and 'pearson'. 'Pearson' assumes normality - check assumptions and make any necessary transformations before using.}

\item{use}{A string specifying either 'pairwise' or 'complete'. 'pairwise' does pairwise deletion of cases and 'complete' selects just the complete cases. Default is 'pairwise'.}

\item{padjust}{A string specifying which adjustment for multiple test should be used. Default is 'BH'. Options are 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', and 'none'.}
}
\value{
A list (or list of lists, if running in groups) containing pca - the pca values for every module and corr - the correlation values for every pca/module against environmental variables.
}
\description{
This function calculates module eigen vectors and then allows for correlation against environmental variables.
}
\details{
A network graph is created for the network(s) that were generated.
}
\examples{
\dontrun{
library(mintJansson)
data(rRNA_data)
mynetwork <- network_calc(omicsData = rRNA_data)
mygraph <- pmartRseq_igraph(netData = mynetwork, coeff=0.6, pval=NULL, qval=0.05)
mymods <- detect_modules(netGraph = mygraph)
myeigen <- mod_eigen(omicsData = rRNA_data, modData=mymods, envVars=c("MBC","MBN","SOC"), method="spearman", use="pairwise", padjust="BH")
myeigen
}

}
\references{
corr.test from psych: Revelle, W. (2017) psych: Procedures for Personality and Psychological Research, Northwestern University, Evanston, Illinois, USA, https://CRAN.R-project.org/package=psych Version = 1.7.8. ;   pca from pcaMethods: Stacklies, W., Redestig, H., Scholz, M., Walther, D. and Selbig, J.  pcaMethods -- a Bioconductor package providing PCA methods for incomplete data. Bioinformatics, 2007, 23, 1164-1167
}
\author{
Allison Thompson
}
