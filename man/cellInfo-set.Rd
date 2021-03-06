% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PharmacoSetClass.R
\name{cellInfo<-}
\alias{cellInfo<-}
\title{cellInfo<- Generic}
\usage{
cellInfo(object) <- value
}
\arguments{
\item{object}{The \code{PharmacoSet} to replace cell info in}

\item{value}{A \code{data.frame} with the new cell annotations}
}
\value{
Updated \code{PharmacoSet}
}
\description{
Generic for cellInfo replace method
}
\examples{
data(CCLEsmall)
cellInfo(CCLEsmall) <- cellInfo(CCLEsmall)

}

