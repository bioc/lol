\name{plotGW}
\alias{plotGW}
\title{
Plot genome-wide data along the genome 
}
\description{
Plot different measurements across the genome such as copy number amplifications and deletions.
}
\usage{
plotGW(data, pos, marks=NULL, fileType='png', file='plotGW', width=1000, height=500,  autoscale=FALSE, col=c('lightblue', 'lightgreen', 'darkblue', 'darkgreen'), legend=1:10, ylab='', pch=19, cex.axis=1.2 ,cex.lab=1.2, cex=.5, legend.pos='bottomright', mtext=NULL, mtext.side=2, mtext.at=NULL, mtext.line=3, ...)}
\arguments{
  \item{data}{
data matrix to plot, each column is plotted individually across the genome
}
  \item{pos}{
the chromosome locations for the data, can be a matrix or data frame with a column named chromosome_name, or a numeric vector
}
  \item{marks}{
if there is specific marks to plot on the baselne, eg. to indicate where are the SNPs, should be a vector of numbers indicating where the marks is relative to the input data matrix
}
  \item{fileType}{
either png or pdf file type
}
\item{file}{file name}
  \item{width}{
width of the plot}
  \item{height}{
height of the plot}
  \item{autoscale}{
should the columns of data be scaled?}
  \item{col}{
colors for each of the data columns to be plotted, should be no shorter than the number of columns in 'data'
}
  \item{legend}{legend text in the legend box 
}
\item{ylab}{parameter for par, default to ''}
\item{pch}{parameter for par, default to 19}
\item{cex.axis}{parameter for par, default to 1.2}
\item{cex.lab}{parameter for par, default to 1.2}
\item{cex}{parameter for par, default to 0.5}
\item{legend.pos}{parameter for legend, default to 'bottomright'}
\item{mtext}{parameter for mtext, default to NULL}
\item{mtext.side}{parameter for mtext, default to 2}
\item{mtext.at}{parameter for mtext, default to 2}
\item{mtext.line}{parameter for mtext, default to 3}
  \item{\dots}{
Other parameters to pass to plot() or legend()
}
}
\details{
This function requires as input data a vector or a matrix with different variables in columns, and a position matrix of chromosome name and start position. The number of rows in the position matrix should be the same as the length of the data vector or the number of rows of the data matrix. The function plots the data according to the position across the genome, providing a genome-wide description.
}
\value{
Write an image file to disk, either in png or pdf format.
}
\author{
Yinyin Yuan
}

\seealso{
lasso.cv
}
\examples{
data(chin07)
gain <- rowSums(chin07$cn >= .2)
loss <- -rowSums(chin07$cn <= -.2)
plotGW(data=cbind(gain, loss), pos=attr(chin07$cn, 'chrome'), legend=c('gain', 'loss'))
}
