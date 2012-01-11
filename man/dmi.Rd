\name{dmi}
\alias{dmi}
\alias{dminjk}
\alias{dmi.pw}
\alias{dminjk.pw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate BCMI for categorical (discrete) data
%%  ~~function to do ... ~~
}
\description{ This function calculates MI and BCMI between a set of discrete variables
    held as columns in a matrix. It also performs jackknife bias correction and
    provides a z-score for the hypothesis of no association. Also included are
    the *.pw functions that calculate MI between two vectors only.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dmi(dmat)
dminjk(dmat)
dmi.pw(disc1, disc2)
dminjk.pw(disc1, disc2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmat}{ The data matrix. Each row is an observation and each column is a
      variable of interest. Should contain categorical data, all types of data
       will be coerced via factors to integers.
}
    \item{disc1}{ A vector for the pairwise version}
    \item{disc2}{ A vector for the pairwise version}
}
\details{ The results of dmi() are in many ways similar to a correlation matrix,
    with each row and column index corresponding to a given variable.
    dminjk() and dminjk.pw() just returns the MI values without performing the
    jackknife.
}
\value{ Returns a list of 3 matrices each of size \code{ncol(dmat)} by
    \code{ncol(dmat)}
    \item{mi}{ The raw MI estimates.}
    \item{bcmi}{ Jackknife bias corrected MI estimates (BCMI). These are each MI value
    minus the corresponding jackknife estimate of bias.}
    \item{zvalues}{ Z-scores for each hypothesis that the corresponding 
    bcmi value is zero. These have poor statistical properties but can be useful
    as a rough measure of the strength of association.}
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
data(cars)

# Discretise the data first
d <- cut(cars$dist, breaks = 10)
s <- cut(cars$speed, breaks = 10)

# Discrete MI values
dmi.pw(s, d)

# For comparison, analysed as continuous data:
cmi.pw(cars$dist, cars$speed)

# Exploring a group of categorical variables
dat <- mtcars[, c("cyl","vs","am","gear","carb")]
discresults <- dmi(dat)
discresults

# Plot the relative magnitude of the BCMI values
diag(discresults$bcmi) <- NA
mp(discresults$bcmi)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
