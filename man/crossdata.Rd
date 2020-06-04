\name{crossdata}
\docType{data}
\alias{crossdata}
\title{The simulated dataset with cumulative incidence function curves of event of interest crossed}
\description{
A dataset constructed by simulation with 200 observations on the following 3 variables.
\describe{
\item{time}{The followed up time of patients.}
\item{status}{The status indicator, 0= censored, 1= event of interest, 2= competing risks.}
\item{group}{The group indicator for comparison, 0= control group, 1= treatment group.}
}
}
\usage{data(crossdata)}
\format{A data frame.}
\keyword{datasets}
