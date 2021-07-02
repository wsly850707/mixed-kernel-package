\name{CauchyP}
\alias{CauchyP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cauchy p-value combination
%%  ~~function to do ... ~~
}
\description{
CauchyP is used to combine the p-values obtained under different w.
}
\usage{
CauchyP(p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
p-values matrix obtained from HDI.
%%     ~~Describe \code{p} here~~
}
}
\details{
P-value obtained from HDI according to different w need to be combined for further inference.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{1-pcauchy(cct)}{
combined p-value vector
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
}
\references{
Liu Y, Xie J. Cauchy Combination Test: A Powerful Test With Analytic p-Value Calculation Under Arbitrary Dependency Structures. J. Am. Stat. Assoc. 2020; 115:393-402.
}
\author{
Yan Liu
}

\keyword{ p-value combination }% use one of  RShowDoc("KEYWORDS")
\keyword{ Cauchy }% __ONLY ONE__ keyword per line