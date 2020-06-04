\name{ABC.ts}
\alias{ABC.ts}
\title{The two-stage test of the Li's test and area between the CIF curves test}
\description{
Statistical inference of two-stage test. Stage I is the Li's test, and stage II is area between the CIF curves(ABC) test.
}
\usage{
ABC.ts(time, status, group, nboot=1000, alpha=0.05, seed=12345)
}
\arguments{
\item{time}{The followed up time for testing data.}
\item{status}{The status indicator, should be coded as 0= censored, 1= event of interest, 2= all other competing events.}
\item{group}{The group indicator for comparison, and the elements of this vector should take either 0 or 1. Normally, 0= control group, 1= treatment group.}
\item{nboot}{The times of bootstrap resamplings, with default as nboot=1000.}
\item{alpha}{The overall significance level, with default as alpha=0.05.}
\item{seed}{The seed number, with default seed=12345.}
}
\value{
\item{method}{Three test results are involved, the Li's test, ABC test, and the two-stage test.}
\item{statistic}{The statistic of the Li's test, ABC test, and the two-stage test.}
\item{Pvalue}{The P value of the Li's test, ABC test, and the two-stage test.}
}
\references{
[1] Li JN, Rademacher JL, Zhang MJ. Weighted comparison of two cumulative incidence functions with R-CIFsmry package. Computer Methods and Programs in Biomedicine, 2014, 116(3): 205-214.

[2] Lyu J, Chen J, Hou Y, Chen Z. Comparison of two treatments in the presence of competing risks. Pharmaceutical Statistics, 2020. DOI: 10.1002/pst.2028.
}
\examples{
#get dataset from package
data(crossdata)
#just for an example, we set resampling times at 10
#two-stage test
ABC.ts(crossdata$time, crossdata$status, crossdata$group, alpha=0.05, nboot=10)
}
