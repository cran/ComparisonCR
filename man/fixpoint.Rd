\name{fixpoint}
\alias{fixpoint}
\title{Statistical inference methods for cumulative incidence functions at a fixed time point}
\description{
A function used to produce the results of three type of statistical inference methods for testing cumulative incidence functions at fixed time point.
}
\usage{
fixpoint(time, status, group, timepoint, type=1)
}
\arguments{
\item{time}{The followed up time for testing data.}
\item{status}{The status indicator, should be coded as 0= censored, 1= event of interest, 2= all other competing events.}
\item{group}{The group indicator for comparison, and the elements of this vector should take either 0 or 1. Normally, 0= control group, 1= treatment group.}
\item{timepoint}{The fixed time point for testing.}
\item{type}{Four type of statistical inference methods. "1" is related to Gaynor variance, "2" is related to Aalen variance. The default is set to type=1.}
}
\value{
\item{method}{Containing 5 transformation of statistical inference methods: line, log, cloglog, arcsin-square, and logist. See more details in references.}
\item{est0}{The estimation of survival rates at the fixed timepoint for control group.}
\item{var0}{The estimation of variances at the fixed timepoint for control group.}
\item{est1}{The estimation of survival rates at the fixed timepoint for treatment group.}
\item{var1}{The estimation of variances at the fixed timepoint for treatment group.}
\item{statistic}{The statistics of corresponding methods.}
\item{Pvalue}{The test P value of corresponding methods.}
}
\references{
Chen J, Hou Y, Chen Z. Statistical inference methods for cumulative incidence function curves at a fixed point in time. Communications in Statistics - Simulation and Computation, 2020, 49(1): 79-94.
}
\examples{
#get dataset from package
data(crossdata)
#
#if there exist differences at 1 year or 3 years
  #Gaynor
fixpoint(crossdata$time, crossdata$status, crossdata$group, timepoint=1, type=1)
  #Aalen
fixpoint(crossdata$time, crossdata$status, crossdata$group, timepoint=3, type=2)
}
