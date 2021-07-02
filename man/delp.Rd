\name{delp}
\alias{delp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
p-value arrangement
%%  ~~function to do ... ~~
}
\description{
delp is used to organize the p-values which are obtained from multiple testing adjustment.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
delp(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  p-values matrix obtained from multiple testing adjustment.
%%     ~~Describe \code{x} here~~
}
}
\details{
If w p-values are all equal to 1 in a variable, then they are not combined by Cauchy. Their combination is 1. If w p-values are all less than 1, then they are combined by Cauchy. If part of p-values are less then 1 and the others are equal to 1, then we just use Cauchy to combine the ones less than 1.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
out
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Yan Liu
}

% R documentation directory.
\keyword{ p-values arrangement }% use one of  RShowDoc("KEYWORDS")
\keyword{ p-values adjustment }% __ONLY ONE__ keyword per line
