\name{matrixLasso}
\alias{matrixLasso}
\title{
A wrapper function for matrix-to-matrix Lasso regressions
}
\description{
This function wraps up different types of lasso optimizers and perform multiple, independent lasso inference on matrix responses. If the dimensionality of the input is small, the function converts the matrix of input response into a vector and solves the problem with one lasso inference. Otherwise, lasso regression is performed independently for each variables in the response matrix. 
}
\usage{
matrixLasso(y, x=NULL, method='cv', nameControl=FALSE, standardize=FALSE, track=0, lambda1=NULL, nFold=10, ...)
}
\arguments{
  \item{y}{
Input response matrix, typically expression data with genes/variables in columns and samples/measurements in rows. Or when input x is NULL, y should be an object of two lists: y: expression data and x: copy number data
}
  \item{x}{
Input predictor matrix, typically copy number data, genes/predictors in columns and samples/measurements in rows. Can be missing if the data is input to y.
}
  \item{method}{
Which optimization method to use for lasso inference, such as 'cv', 'stability', 'simultaneous', and 'multiSplit'.
}
  \item{nameControl}{
If the same item appears in both responses and predictors, the regression should remove the one same as the response from the predictors. This happens when for example a single data type is use for inferring gene network from expression data. Enable nameControl in this case. 
}
  \item{standardize}{
Option to standardize the data, default to TRUE
}
  \item{track}{
Option to display progress, default to 0, 1 gives a brief summary of each fit, and 2 gives the full detail.
}
  \item{lambda1}{
The minimum lambda to use, default to NULL for which the program will select it automatically
}
  \item{nFold}{
Number of folds for cross-validation, default to 10
}
  \item{\dots}{
}
}
\value{
\item{coefMat}{A coefficient matrix, rows are responses and columns are predictors}
\item{fit}{If only a single regression is used for matrix lasso, the fit return.}
\item{resMat}{A residual matrix, each row is the residuals of a response.}
}
\author{
Yinyin Yuan
}
\seealso{
lasso
}
\examples{
data(chin07)
data <- list(y=t(chin07$ge), x=t(chin07$cn))
res <- matrixLasso(data, method='cv', nFold=5)
res
}