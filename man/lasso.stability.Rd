\name{lasso.stability}
\alias{lasso.stability}
\title{
Stability and randomised lasso
}
\description{
point-wise controled lasso stability selection
}
\usage{
lasso.stability(y, x=NULL, alpha=.5, subsampling=.5, nSubsampling=200, model='linear', pi_th=.6,  alpha.fwer=1, lambda1=NULL, steps=10, track=FALSE,  standardize=FALSE,  ...)
}
\arguments{
  \item{y}{A vector of gene expression of a probe, or a list object if x is NULL. In the latter case y should a list of two components y and x, y is a vector of expression and x is a matrix containing copy number variables}
  \item{x}{Either a matrix containing CN variables or NULL}
 \item{alpha}{weakness parameter: control the shrinkage of regulators, if alpha = 1 then no randomisation, if NULL then a randomly generated vector is used}
  \item{subsampling}{fraction  of samples to use in the sampling process, default to 0.5}
  \item{nSubsampling}{The number of subsampling to do, default to 200
}
  \item{model}{which model to use, one of "cox", "logistic", "linear",
         or "poisson". Default to 'linear'}
  \item{pi_th}{
The threshold of the stability probablity for selecting a regulator. It is to determine whether a coefficient is non-zero based on the frequency it is subsampled to be non-zero, default to 0.6
}
  \item{alpha.fwer}{Parameter to control for the FWER, choosing alpha.fwer and alpha control the E(V), V being the number of noise variables, eg. when alpha=0.9, alpha.fwer = 1 control the E(V)<=1}
  \item{lambda1}{
minimum lambda to use}
  \item{steps}{
parameter to be passed on to penalized}
  \item{track}{
track the progress, 0 none tracking, 1 minimum amount of information and 2 full information}
  \item{standardize}{
standardize the data or not?
}
  \item{\dots}{
}
}
\details{
The function first selects lambda that approximately give maximum sqrt(.8*p) predictors, while p is the number of total predictors. Then it runs lasso a number of times keeping lambda fixed. These runs are randomised with scaled predictors and subsamples. At the end, the non-zero coefficients are determined by their frequencies of selections. 
}
\value{
A list object of class 'lol', consisting of:
\item{beta}{coefficients}
\item{beta.bin}{binary beta vector as thresholded by pi_th}
\item{mat}{the sampling matrix, each column is the result of one sampling}
\item{residuals}{residuals of regression model}
}
\references{
N. Meinshausen and P. Buehlmann (2010), Stability Selection (with discussion), Journal of the Royal Statistical Society, Series B, 72, 417-473.

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
res <- lasso.stability(data, nSubsampling=50)
res
}