% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemFunctions.R
\name{getMolecularFormula}
\alias{getMolecularFormula}
\title{Return a molecular formula}
\usage{
getMolecularFormula(form = "", output = "table", order = "hill")
}
\arguments{
\item{form}{Input formula}

\item{output}{output type - either table (columns = Element and Count) or formula (string)}

\item{order}{atom ordering - either alpha or hill notation}
}
\value{
a tibble or string corresponding to the molecular formula
}
\description{
Return a standardized molecular formula
}
