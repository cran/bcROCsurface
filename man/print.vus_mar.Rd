% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcvusMAR.R
\name{print.vus_mar}
\alias{print.vus_mar}
\title{Print summary results of VUS}
\usage{
\method{print}{vus_mar}(x, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
\item{x}{an object of class "vus_mar", a result of a call to \code{\link{vus_mar}}.}

\item{digits}{minimal number of significant digits, see \code{\link{print.default}}.}

\item{...}{further arguments passed to \code{\link{print}} method.}
}
\description{
\code{print.vus_mar} prints the results for the output of function \code{\link{vus_mar}}.
}
\details{
\code{print.vus_mar} shows a nice format of the summary table for the VUS estimate results. Some information on the diagnostic test, the fitted values of VUS, and confidence intervals are shown.
}
\seealso{
\code{\link{vus_mar}}
}
