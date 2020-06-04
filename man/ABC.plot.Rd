\name{ABC.plot}
\alias{ABC.plot}
\title{Plot for area between the CIF curves}
\description{
A function used to plot the area between the CIF curves.
}
\usage{
ABC.plot(time, status, group, tau=NULL, max.x=NULL, max.y=NULL, col=c(1,1,8),lwd=c(3,3)
         , lty=c(1,2), lab.x="Time", lab.y="CIF", cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
}
\arguments{
\item{time}{The followed up time for testing data.}
\item{status}{The status indicator, should be coded as 0= censored, 1= event of interest, 2= all other competing events.}
\item{group}{The group indicator for comparison, and the elements of this vector should take either 0 or 1. Normally, 0= control group, 1= treatment group.}
\item{tau}{The truncation time point for shading lines, which needs to be smaller than or equal to the maximum of time. When tau=NULL, the default value is used. See details in reference.}
\item{max.x}{The maximum of x-axis. When the default NULL is used, the max.x is set to the max velue of time.}
\item{max.y}{The maximum of y-axis. When the default NULL is used, the max.y is set to 1.}
\item{col}{A vector of integers specifying colors for control group, treatment group, and shading lines, respectively. The default value is c("black","black","gray90").}
\item{lwd}{A vector of numeric values of line widths for control group, treatment group, and shading lines, respectively.. The default value is c(3,3).}
\item{lty}{A vector of integers specifying line types for control group, treatment group, and shading lines, respectively. The default value is c(1,2).}
\item{lab.x}{Label given to the x-axis.}
\item{lab.y}{Label given to the y-axis.}
\item{cex.axis}{The magnification to be used for axis annotation relative to the current setting of cex.}
\item{cex.lab}{The magnification to be used for x and y labels relative to the current setting of cex.}
\item{cex.main}{The magnification to be used for main titles relative to the current setting of cex.}
}
\value{None}
\references{
Lyu J, Chen J, Hou Y, Chen Z. Comparison of two treatments in the presence of competing risks. Pharmaceutical Statistics, 2020. DOI: 10.1002/pst.2028.
}
\note{
ABC.plot is based on plot function and polygon function.
}
\seealso{
\code{\link{ABC}}
}
\examples{
#get dataset from package
data(crossdata)
#
#get plot for ABC
ABC.plot(crossdata$time, crossdata$status, crossdata$group, col=c(1,1,"gray90"))
legend(0,1.02,c("group0","group1"), col=c(1,1), lty=c(1,2), lwd=c(3,3), cex=1.2, bty="n")
}
