% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeIC50.R
\name{computeIC50}
\alias{computeIC50}
\title{Fits dose-response curves to data given by the user
and returns the IC50 of the fitted curve.}
\usage{
computeIC50(conc, viability, trunc = TRUE, verbose = FALSE)
}
\arguments{
\item{conc}{[vector] is a vector of drug concentrations.}

\item{viability}{[vector] is a vector whose entries are the viability values observed in the presence of the
drug concentrations whose logarithms are in the corresponding entries of the log_conc, expressed as percentages
of viability in the absence of any drug.}

\item{trunc}{[logical], if true, causes viability data to be truncated to lie between 0 and 1 before
curve-fitting is performed.}

\item{verbose}{[logical] should diagnostic messages be printed? (default=FALSE)}
}
\value{
An estimate of the IC50 for the concentrations and viabilities provided
}
\description{
Fits dose-response curves to data given by the user
and returns the IC50 of the fitted curve.
}
\examples{
dose <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8") 
viability <- c("108.67","111","102.16","100.27","90","87","74","57")
computeIC50(dose, viability)

}

