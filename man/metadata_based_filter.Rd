% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadata_based_filter.R
\name{metadata_based_filter}
\alias{metadata_based_filter}
\title{Metadata filter object}
\usage{
metadata_based_filter(omicsData, criteria)
}
\arguments{
\item{omicsData}{An object of the classes "seqData"}

\item{criteria}{Specify which omicsData$e_meta column name to filter on}
}
\value{
An object of class metaFilter (also a data.frame) that contains the
  sample identifier, the values in the criteria column, and the corresponding
  sum abundances.
}
\description{
This function returns a metaFilter object, performing all of the
calculations needed to filter the data based off a specified
e_meta level
}
\examples{

}
\author{
Allison Thompson and Sarah Reehl
}
