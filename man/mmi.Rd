\name{mmi}
\alias{mmi}
\alias{mminjk}
\alias{mmi.pw}
\alias{mminjk.pw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate mixed-pair BCMI between a set of continuous variables and a set
    of discrete variables.
%%  ~~function to do ... ~~
}
\description{This function calculates MI and BCMI between a set of continuous variables
    and a set of discrete variables (variables in columns). It also performs
    jackknife bias correction and provides a z-score for the hypothesis of no
    association. Also included are the *.pw functions that calculate MI between
    two vectors only. The *njk functions do not perform the jackknife and are
    therefore faster. 
    %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mmi(cts, disc, level = 3L, na.rm = FALSE, h, ...)
mminjk(cts, disc, level = 3L, na.rm = FALSE, h, ...)
mmi.pw(cts, disc, h, ...) 
mminjk.pw(cts, disc, h, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cts}{ The data matrix. Each row is an observation and each column is a
      variable of interest. Should be numerical data. (For the pairwise functions this
	should be a vector.)
%%     ~~Describe \code{cts} here~~
}
    \item{disc}{ Matrix of discrete data, each row is an observation and each
        column is a variable. Will be coerced to integers. (For the pairwise functions this
	should be a vector.)}
%%
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
}
\details{mminjk() and mminjk.pw() return just the MI values without performing the
    jackknife. mmi.pw() and mminjk.pw() only require one bandwidth for the
    continuous variable. The number of processor cores used can be changed by
    setting the environment variable "OMP_NUM_THREADS" \emph{before} starting R.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns a list of 3 matrices each of size ncol(cts) by
    ncol(disc). Each row index represents a continuous variable and each column
    index a discrete variable. 
    \item{mi}{The raw MI estimates.}
    \item{bcmi}{Jackknife bias corrected MI estimates (BCMI). These are each MI value
    minus the corresponding jackknife estimate of bias.}
    \item{zvalues}{z-scores for each hypothesis that the corresponding 
    bcmi value is zero. These have poor statistical properties but can be useful
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
# A dataset with discrete and continuous variables

cts <- state.x77
disc <- data.frame(state.division,state.region)
summary(cts)
table(disc)
m1 <- mmi(cts, disc)
lapply(m1, round, 2)
# Division gives more information about the continuous variables than region.

# Here is one where both division and region show a strong association:
boxplot(cts[,6] ~ disc[,1])
boxplot(cts[,6] ~ disc[,2])

# In this case the states need to be divided into regions before a clear
# association can be seen:
boxplot(cts[,1] ~ disc[,1])
boxplot(cts[,1] ~ disc[,2])

# Look at associations within the continuous variables:
pairs(cts, col = state.region)
c1 <- cmi(cts)
lapply(c1, round, 2)

##################################################
# A pairwise comparison

# Note that the ANOVA homoskedasticity assumption is not satisfied here.
boxplot(InsectSprays[,1] ~ InsectSprays[,2])
mmi.pw(InsectSprays[,1], InsectSprays[,2])

##################################################
# Another pairwise comparison

boxplot(morley[,3] ~ morley[,1])
m2 <- mmi.pw(morley[,3], morley[,1])
m2

##################################################
# See the vignette for large-scale examples.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
