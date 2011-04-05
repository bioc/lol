\name{lol-package}
\alias{lol-package}
\alias{lol}
\docType{package}
\title{
Lots of Lasso}
\description{
Various optimization methods for Lasso inference with matrix wrapper. 
}
\details{
\tabular{ll}{
Package: \tab lol\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2011-04-02\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

}
\author{
Yinyin Yuan
Maintainer: Yinyin Yuan <yy341@cam.ac.uk>
}
\references{
Goeman, J. J. (2009), L1 penalized estimation in the cox proportional hazards model. Biometrical Journal.
N. Meinshausen and P. Buehlmann (2010), Stability Selection (with discussion), Journal of the Royal Statistical Society, Series B, 72, 417-473.
}
\keyword{ package }
\seealso{
lasso, matrixLasso
}
\examples{
data(chin07)
data <- list(y=t(chin07$ge), x=t(chin07$cn))
res <- matrixLasso(data, method='cv', nFold=5)
res
}
