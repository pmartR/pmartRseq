% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm_funcs.R
\name{Log_Norm}
\alias{Log_Norm}
\title{Log transformation of count data}
\usage{
Log_Norm(e_data, edata_id)
}
\arguments{
\item{e_data}{a \eqn{p \times n} data.frame of count data, where \eqn{p} is the number of features and \eqn{n} is the number of samples. Each row corresponds to data for a feature, with the first column giving the feature name.}

\item{edata_id}{character string indicating the name of the feature identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.}
}
\value{
List containing 3 elements: norm_data is a data.frame with same structure as e_data that contains the Log2-normalized data, location_param is NULL, scale_param is a numeric vector of 1's, for later use.
}
\description{
The method normalizes count data by log2 transforming all of the counts
}
\details{
Count data is normalized by a log2 transformation
}
\examples{
library(mintJansson)
data(cDNA_hiseq_data)
cDNA_log <- Log_Norm(e_data = cDNA_hiseq_data$e_data, edata_id = attr(cDNA_hiseq_data, "cnames")$edata_cname)

}
\author{
Allison Thompson
}
