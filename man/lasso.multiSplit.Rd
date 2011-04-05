\name{lasso.multiSplit}
\alias{lasso.multiSplit}
\title{
Multi-split lasso 
}
\description{
Multi-split lasso as described in Meinshausen 2009
}
\usage{
lasso.multiSplit(y, x=NULL, lambda1=NULL, nSubsampling=200, model='linear', alpha=0.05, gamma.min=0.05, gamma.max=0.95, track=FALSE, ...)
}
\arguments{
  \item{y}{A vector of gene expression of a probe, or a list object if x is NULL. In the latter case y should a list of two components y and x, y is a vector of expression and x is a matrix containing copy number variables}
  \item{x}{Either a matrix containing CN variables or NULL}
  \item{nSubsampling}{number of splits, default to 200}
  \item{model}{which model to use, one of "cox", "logistic", "linear",
         or "poisson". Default to 'linear'}
  \item{alpha}{specify significant level to determine the non-zero coefficients in the range of 0 and 1, default to 0.05}
  \item{gamma.min}{
the lower bound of gamma
}
  \item{gamma.max}{
the higher bound of gamma
}
  \item{lambda1}{
minimum lambda to be used, if known
}
  \item{track}{
track progress}
  \item{\dots}{
other parameters to be passed to lass.cv}
}
\details{
This function performs the multi-split lasso as proposed by Meinshausen et al. 2009. The samples are first randomly split into two disjoint sets, one of which is used to find non-zero coefficients with a regular lasso regression, then these non-zero coefficients are fitted to another sample set with OLS. The resulting p-values after multiple runs can then be aggregated using quantiles. 
}
\value{
A list object of class 'lol', consisting of:
\item{beta}{coefficients}
\item{mat}{the Q_gamma matrix as described in the paper}
\item{residuals}{residuals, here is only the input y}
\item{pmat}{the adjusted p matrix as described in the paper}
}
\references{
Nicolai Meinshausen, Lukas Meier and Peter Buehlmann (2009),
P-values for high-dimensional regression. 
Journal of the American Statistical Association, 104, 1671-1681.
}
\author{
Yinyin Yuan
}
\seealso{
lasso
}
\examples{
data(chin07)
data <- list(y=chin07$ge[1,], x=t(chin07$cn))
res <- lasso.multiSplit(data, nSubsampling=50)
res
}