\name{lmMatrixFit}
\alias{lmMatrixFit}
\title{
Multiple lm fit for penalized regressions
}
\description{
Refit the regressions given matrices of responses, predictors, and the coefficients/interactions matrix. This is typically used after the lasso, since the coefficients were shrinked.
}
\usage{
lmMatrixFit(y, x = NULL, mat, th = NULL)
}
\arguments{
  \item{y}{
Input response matrix, typically expression data with genes/variables in columns and samples/measurements in rows.  Or when input x is NULL, y should be an object of two lists: y: expression data and x: copy number data
}
  \item{x}{
Input predictor matrix, typically copy number data, genes/predictors in columns and samples/measurements in rows. Can be NULL
}
  \item{mat}{
Coefficient matrix, number of columns is the number of predictors (y) and number of rows is the number of responses (x)
}
  \item{th}{
The threshold to use in order to determine which coefficients are non-zero, so the corresponding predictors are used
}
}
\value{
\item{coefMat}{A coefficient matrix, rows are responses and columns are predictors}
\item{resMat}{A residual matrix, each row is the residuals of a response.}
\item{pvalMat}{Matrix of p-values for each coefficients}
}
\author{
Yinyin Yuan
}
\seealso{
lm, matrixLasso
}
\examples{
data(chin07)
data <- list(y=t(chin07$ge), x=t(chin07$cn))
res <- matrixLasso(data, method='cv', nFold=5)
res
res.lm <- lmMatrixFit(y=data, mat=abs(res$coefMat), th=0.01)
res.lm
}