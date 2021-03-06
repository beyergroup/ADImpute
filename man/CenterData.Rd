% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NetworkFunctions.R
\name{CenterData}
\alias{CenterData}
\title{Data centering}
\usage{
CenterData(data)
}
\arguments{
\item{data}{matrix of gene expression to be centered row-wise (genes as rows
and samples as columns)}
}
\value{
list; row-wise centers and centered data
}
\description{
\code{CenterData} centers expression of each gene at 0
}
