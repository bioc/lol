\name{lasso.simultaneous}
\alias{lasso.simultaneous}
\title{
Simultaneous lasso
}
\description{
The function performs lasso with multiple random sample splits, selecting coefficients that are simultaneously non-zero in both subsets of samples.
}
\usage{
lasso.simultaneous(y, x=NULL, model='linear', nSubsampling=200, alpha=.5, lambda1=NULL, track=FALSE, ...)
}
\arguments{
  \item{y}{A vector of gene expression of a probe, or a list object if x is NULL. In the latter case y should a list of two components y and x, y is a vector of expression and x is a matrix containing copy number variables}
  \item{x}{Either a matrix containing CN variables or NULL}
  \item{model}{which model to use, one of "cox", "logistic", "linear",
         or "poisson". Default to 'linear'}
  \item{nSubsampling}{
The number of random permutations, both on sample spliting and on variable scaling, default to 200.
}
 \item{alpha}{weakness parameter: control the shrinkage of regulators. The lower alpha is, the bigger the vanishing effect on small coefficients.}
\item{lambda1}{minimum lambda, default to NULL}
\item{track}{logical value, whether to track the progress}
  \item{\dots}{
Other parameters to be passed to the penalized function
}
}
\details{
In each run the function splits samples randomly to two equal sets, run lasso on both sets, then select those coefficients that are simultaneously non-zero across two sets. Finally the results across many runs are summarized as the frequency of selected predictors - the higher the frequency the more confidence that the corresponding predictors are significant.
}
\value{
A list object of class 'lol', consisting of:
\item{beta}{Coefficient vector}
\item{n}{Number of actual subsampling, should be equal or smaller than nSubsampling in case of failing.}
\item{mat}{result matrix of the subsampling}
}
\references{
N. Meinshausen and P. Buehlmann (2010), Stability Selection (with discussion), Journal of the Royal Statistical Society, Series B, 72, 417-473.}
\author{
Yinyin Yuan
}
\seealso{
lasso
}
\examples{
data(chin07)
data <- list(y=chin07$ge[1,], x=t(chin07$cn))
res <- lasso.simultaneous(data, nSubsampling=50)
res
}