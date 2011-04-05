\name{lasso}
\alias{lasso}
\title{
lasso
}
\description{
Lasso penalized linear regression with different optimizers
}
\usage{
lasso(y, ...)
}
\arguments{
  \item{y}{A list object of one of the four classes: 'cv', 'stability',  'multiSplit', and 'simultaneous'. If x is NULL then y should a list of two components y and x, y is a vector of expression and x is a matrix containing copy number variables}
  \item{...}{other parameters}
}
\details{
The function contains various optimization methods for Lasso inference, such as cross-validation, randomised lasso, simultaneous lasso etc. It is specifically designed for multicollinear predictor variables. 
}
\value{
Varied depending on the optimizer used. Generally it contains
\item{beta}{coefficients}
\item{residuals}{residuals of regression model}
\item{fit}{the corresponding fit of regression}
}
\references{
Goeman, J. J. (2009), L1 penalized estimation in the cox proportional hazards model, Biometrical Journal.
N. Meinshausen and P. Buehlmann (2010), Stability Selection (with discussion), Journal of the Royal Statistical Society, Series B, 72, 417-473.
Nicolai Meinshausen, Lukas Meier and Peter Buehlmann (2009),
P-values for high-dimensional regression. 
Journal of the American Statistical Association, 104, 1671-1681.

}
\author{
Yinyin Yuan
}
\seealso{matrixLasso
}
\examples{
data(chin07)
data <- list(y=chin07$ge[1,], x=t(chin07$cn))
class(data) <- 'cv'
res <- lasso(data)
}