\name{chin07}
\alias{chin07}
\docType{data}
\title{
Breast cancer data set of genome-wide copy number merged data and expression of some important genes
}
\description{
A subset of breast cancer data as used in Yuan et al. (to be submitted).
}
\usage{data(chin07)}
\format{
A list object of two named data matrices, cn: DNA copy number, ge: RNA expression. The matrices columns are samples and rows are probes/variables.
}
\details{
Genome-wide copy number data was merged using CGHregions resulting in 339 regions across 106 samples. Expression data are 7 probes mapped to important breast cancer genes such as CCNE2, MYC, etc, also of 106 samples.
}
\references{
Chin SF, Teschendorff AE, Marioni JC, Wang Y, Barbosa-Morais NL, et al. (2007) High-resolution arraycgh and expression profiling identifies a novel genomic subtype of er negative breast cancer. Genome Biology 8: R215+.
Yuan et al. (2011) Discovery and functional annotation of cis- and trans-acting DNA copynumber hotspots in breast cancer, to be submitted.
}
\examples{
data(chin07)
gain <- rowSums(chin07$cn >= .2)
loss <- -rowSums(chin07$cn <= -.2)
plotGW(data=cbind(gain, loss), pos=attr(chin07$cn, 'chrome'), legend=c('gain', 'loss'))
}
\keyword{datasets}
