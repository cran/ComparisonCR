\name{ABC.perm}
\alias{ABC.perm}
\title{The permutation test of area between the CIF curves for competing risk data}
\description{
Statistical inference of area between the CIF curves (ABC) test with a permutation procedure for competing risk data.
}
\usage{
ABC.perm(time, status, group, t0=0, tau=NULL, nperm=1000, seed=12345, bias=FALSE)
}
\arguments{
\item{time}{The followed up time for testing data.}
\item{status}{The status indicator, should be coded as 0= censored, 1= event of interest, 2= all other competing events.}
\item{group}{The group indicator for comparison, and the elements of this vector should take either 0 or 1. Normally, 0= control group, 1= treatment group.}
\item{t0}{The start time point for testing, the default is set as t0=0 for overall test.}
\item{tau}{The truncation time point for testing, which needs to be smaller than or equal to the minimum of the largest observed time in each of the two groups. When tau=NULL, the default value is used. See details in reference.}
\item{nperm}{The times of permutation, with default nperm=1000.}
\item{seed}{The seed number, with default seed=12345.}
\item{bias}{If bias=TRUE, the results will print the mean difference for delta of simulated datasets and original data. The default is bias=FALSE}
}
\value{
\item{t0}{The start time point for testing.}
\item{tau}{The truncation time point for testing.}
\item{delta}{The alsolute difference of two cumulative incidence functions.}
\item{var(delta)}{The variance for delta of simulated datasets.}
\item{bias}{The mean difference for delta of simulated datasets and original data.}
\item{Pvalue}{The P value of this test.}
}
\references{
Lyu J, Chen J, Hou Y, Chen Z. Comparison of two treatments in the presence of competing risks. Pharmaceutical Statistics, 2020. DOI: 10.1002/pst.2028.
}
\examples{
#get dataset from package
data(crossdata)
#just for an example, we set resampling times at 10
#overall test
ABC.perm(crossdata$time, crossdata$status, crossdata$group, nperm=10)
#arbitary test for detecting difference after 2 years
ABC.perm(crossdata$time, crossdata$status, crossdata$group, t0=2, tau=NULL, nperm=10)
#arbitary test for detecting difference between 2 years and 4 years
ABC.perm(crossdata$time, crossdata$status, crossdata$group, t0=2, tau=4, nperm=10)
}
