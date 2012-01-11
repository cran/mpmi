\name{cmi}
\alias{cmi}
\alias{cminjk}
\alias{cmi.pw}
\alias{cminjk.pw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate BCMI between a set of continuous variables
%%  ~~function to do ... ~~
}
\description{This function calculates MI and BCMI between a set of continuous variables
    held as columns in a matrix. It also performs jackknife bias correction and
    provides a z-score for the hypothesis of no association. Also included are
    the *.pw functions that calculate MI between two vectors only. The *njk
    functions do not perform the jackknife and are therefore faster. 
    %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cmi(cts, level = 3L, na.rm = FALSE, h, ...)
cminjk(cts, level = 3L, na.rm = FALSE, h, ...)
cmi.pw(v1, v2, h, ...)
cminjk.pw(v1, v2, h, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cts}{ The data matrix. Each row is an observation and each column is a
      variable of interest. Should be numerical data.
%%     ~~Describe \code{cts} here~~
}
  \item{level}{ The number of levels used for plug-in bandwidth estimation (see
      the documentation for the KernSmooth package.)
%%     ~~Describe \code{level} here~~
}
  \item{na.rm}{ Remove missing values if TRUE. This is required for the
      bandwidth calculation.
%%     ~~Describe \code{na.rm} here~~
}
  \item{h}{ A (double) vector of smoothing bandwidths, one for each variable. If
      missing this will be calculated using the dpik() function from the
      KernSmooth package.
%%     ~~Describe \code{h} here~~
}
  \item{\dots}{ Additional options passed to dpik() if necessary.
%%     ~~Describe \code{\dots} here~~
} 
    \item{v1}{ A vector for the pairwise version}
    \item{v2}{ A vector for the pairwise version}
}
\details{ The results of cmi() are in many ways similar to a correlation matrix,
    with each row and column index corresponding to a given variable.
    cminjk() and cminjk.pw() just returns the MI values without performing the
    jackknife. cmi.pw() and cminjk.pw() each only require two bandwidths, one
    for each variable. The number of processor cores used can be changed by
    setting the environment variable "OMP_NUM_THREADS" \emph{before} starting R.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns a list of 3 matrices each of size ncol(cts) by
    ncol(cts)
    \item{mi}{ The raw MI estimates.}
    \item{bcmi}{Jackknife bias corrected MI estimates (BCMI). These are each MI value
    minus the corresponding jackknife estimate of bias.}
    \item{zvalues}{Z-scores for each hypothesis that the corresponding 
    BCMI value is zero. These have poor statistical properties but can be useful
    as a rough measure of the strength of association.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
%% \author{
%% %%  ~~who you are~~
%% }
%% \note{
%% %%  ~~further notes~~
%% }
%% 
%% %% ~Make other sections like Warning with \section{Warning }{....} ~
%% 
%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
##################################################
# The USArrests dataset

# Matrix version
c1 <- cmi(USArrests)
lapply(c1, round, 2)

# Pairwise version
cmi.pw(USArrests[,1], USArrests[,2])

# Without jackknife
c2 <- cminjk(USArrests)
round(c2, 2)
cminjk.pw(USArrests[,1], USArrests[,2])

##################################################
# A look at Anscombe's famous dataset.
par(mfrow = c(2,2))
plot(anscombe$x1, anscombe$y1)
plot(anscombe$x2, anscombe$y2)
plot(anscombe$x3, anscombe$y3)
plot(anscombe$x4, anscombe$y4)

cor(anscombe$x1, anscombe$y1)
cor(anscombe$x2, anscombe$y2)
cor(anscombe$x3, anscombe$y3)
cor(anscombe$x4, anscombe$y4)

cmi.pw(anscombe$x1, anscombe$y1)
cmi.pw(anscombe$x2, anscombe$y2)
cmi.pw(anscombe$x3, anscombe$y3)
# dpik() has some trouble with zero scale estimates on this one:
cmi.pw(anscombe$x4, anscombe$y4, scalest = "stdev")
##################################################

##################################################
# The highly collinear Longley dataset

pairs(longley, main = "longley data")
l1 <- cmi(longley)
lapply(l1, round, 2)

# Here we demonstrate the scale-invariance of MI.
# Note: Scaling can help stabilise estimates when there are
# difficulties with the bandwidth estimation, but is unnecessary
# here.
long2 <- scale(longley)
l2 <- cmi(long2)
lapply(l2, round, 2)

##################################################
# See the vignette for large-scale examples.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
