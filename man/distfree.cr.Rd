\name{distfree.cr}
\alias{distfree.cr}
\title{
distfree.cr
}
\description{
Constructs empirical confidence regions for bivariate data based on the method proposed by Zhiqiu Hu and Rong-cai Yang(2013) <doi:10.1371/journal.pone.0081179.g001>.
}
\usage{
distfree.cr(x, y, alpha = 0.05, alpha.min.diff = 0.5/NROW(x), nknots = 40, 
            xlab = deparse(substitute(x)), ylab = deparse(substitute(y)), 
            col = c("red", "black", "gray"), draw = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric vector, of dimensions \emph{nobs} * 1. If a data frame or a two-column numeric matrix of x and y is supplied here, the second option \code{y} of the function needs to be ignored.}
  \item{y}{numeric vector, of dimensions \emph{nobs} * 1. This option needs to be ignored if users provided both x and y in the first option of the function.}
  \item{alpha}{Significant level. By default \code{alpha} is set to be 0.05.}
  \item{alpha.min.diff}{minimum difference is allowed for calculating confidence region. This option is not suggested for most users. The default value is set to be \code{alpha/10.}}
  \item{nknots}{number of knots that will be used to enclose the confidence region. The default value \code{nknots=40} is recommended for all users.}
  \item{xlab}{define the label of x axis of the plot.}
  \item{ylab}{define the label of y axis of the plot.}
  \item{col}{define colors of the scatter points and lines of the plot. The default setting \code{col=("red", "black", "gray")} are the colors for the lines enclosed the region, the points within the region and the points outside of the region, respectively.}
  \item{draw}{a logical indicator. Users may disable plotting by setting the option to \code{FALSE}}
}
\details{
This function constructs a distribution-free confidence region based on the method proposed by Zhiqiu Hu and Rong-cai Yang.
}
\value{

\item{alpha.realized}{Realized-alpha, which is defined as the proportion of the total points outside the confidence region.}
\item{polygon}{'data.frame' of x,y providing the apexes of the lines.}
\item{polygon.smooth1}{'data.frame' of x,y providing the apexes of the smoothed polygon 1.}
\item{polygon.smooth2}{'data.frame' of x,y providing the apexes of the smoothed polygon 2.}
\item{data}{'data.frame', of dimension \emph{nobs} * 3, the first two columns are input data of x and y values and the third column \code{data$pip} are indicators of whether the points are within (\code{1}) or outside (\code{0}) the confidence region.}
\item{alpha,xlab,ylab,col}{values assigned by users.}
}

\author{
Zhiqiu Hu and Rong-cai Yang
}
\note{
A smooth confidence region can be achieved by setting up a big number for input variable nknots, and this in turn requires large sample sizes. 
}

\examples{
library(distfree.cr)
dat=data.frame(x=c(rnorm(3000), rnorm(3000, mean=1, sd=2.5)), 
               y=c(rnorm(3000), rnorm(3000, mean=1, sd=2.5)))
pt=distfree.cr(dat, draw=TRUE, alpha=0.05)
pt=distfree.cr(x=dat$x, y=dat$y, draw=FALSE)
plot(pt)
}