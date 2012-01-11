\name{mp}
\alias{mp}
\title{Matrix Plot}
\description{
    Plot a matrix of values in the same order that it is stored 
(the usual mathematical way).
}
\usage{
mp(mat, ...)
}
\arguments{
  \item{mat}{A numeric matrix to be plotted}
  \item{...}{Additional arguments to pass to \code{image()}}
}
\details{
The \code{mp()} function is a simple wrapper to
\code{image()} with a few minor changes. The plot is
flipped so that points correspond to their position in the
matrix. Also, the colours are scaled so that red is the
maximum value in the matrix and white is the minimum value.
}
\examples{
# From the vignette:
data(mpmidata)
ctsresult <- cmi(cts)
mp(ctsresult$bcmi)
} 
