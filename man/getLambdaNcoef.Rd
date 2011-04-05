\name{getLambdaNcoef}
\alias{getLambdaNcoef}
\title{
get the lambda value that yield certain number of non-zero coefficients}
\description{
get the lambda value that yield certain number of non-zero coefficients
}
\usage{
getLambdaNcoef(y, x, lambda1, nCoef, track=FALSE, model='linear', standardize=FALSE)
}
\arguments{
  \item{y}{A vector of expressions}
  \item{x}{a matrix of CN variables}
  \item{lambda1}{
minimum lambda to use}
  \item{nCoef}{
the number of coefficients to get}
  \item{track}{
logical value for tracking the progress}
  \item{model}{which model to use, default to 'linear'
}
  \item{standardize}{
standardize the data or not
}
}
\value{
\item{lambda}{The lambda value that gives approximate same number of non-zero coefficients as required}
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
getLambdaNcoef(data$y, data$x, lambda1=.1, nCoef=10, track=TRUE)
}