\name{lasso.cv}
\alias{lasso.cv}
\title{
Cross validation optimizer for lasso
}
\description{
Cross validation lasso. This function optimizes the lasso solution for correlated regulators by an algorithm. this algorithm chooses the minimum lambda since the penalized package by default use 0 for the minimum, which sometimes take a long time to compute
}
\usage{
lasso.cv(y, x=NULL, lambda1=NULL, model='linear', steps=15, minsteps=5, log=TRUE, track=FALSE, standardize= FALSE, unpenalized=~0, nFold=10, nMaxiter = Inf, ...)
}
\arguments{
  \item{y}{A vector of gene expression of a probe, or a list object if x is NULL. In the latter case y should a list of two components y and x, y is a vector of expression and x is a matrix containing copy number variables}
  \item{x}{Either a matrix containing CN variables or NULL}
  \item{lambda1}{
minimum lambda to use}
  \item{model}{which model to use, one of "cox", "logistic", "linear",
         or "poisson". Default to 'linear'}
  \item{steps}{
parameter to be passed to penalized
}
  \item{minsteps}{
parameter to be passed to penalized
}
  \item{log}{
parameter to be passed to penalized
}
  \item{track}{
parameter to be passed to penalized
}
  \item{standardize}{
parameter to be passed to penalized
}
  \item{unpenalized}{
parameter to be passed to penalized
}
  \item{nFold}{
parameter to be passed to penalized
}
  \item{nMaxiter}{
parameter to be passed to penalized
}
\item{...}{other parameter to be passed to penalized
}
}
\value{
A list object of class 'lol', consisting of:
\item{fit}{The final sparse regression fit}
\item{beta}{the coefficients, non-zero ones are significant}
\item{lambda}{the penalty parameter lambda used}
\item{residuals}{regression residuals}
\item{conv}{logical value indicating whether the optimization has converged}
}
\references{
Goeman, J. J. (2009), L1 penalized estimation in the cox proportional hazards model, Biometrical Journal.
}
\author{
Yinyin Yuan
}
\seealso{lasso
}
\examples{
data(chin07)
data <- list(y=chin07$ge[1,], x=t(chin07$cn), nFold=5)
res <- lasso.cv(data)
res
}