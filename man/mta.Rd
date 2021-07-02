\name{mta}
\alias{mta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
multiple testing adjustment
%%  ~~function to do ... ~~
}
\description{
mta is used to adjust the p-values obtained from multiple tests.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mta(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  main and interaction effects matrix that is used to HDI.
%%     ~~Describe \code{x} here~~
}
}

\value{
 \item{sum.fx}{
effective number of independent tests
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
}
\references{
 Li J, Ji L. Adjusting multiple testing in multilocus analyses using the eigenvalues of a correlation matrix. Heredity (Edinb). 2005; 95:221-227.
%% ~put references to the literature/web site here ~
}
\author{
Yan Liu
}

\keyword{ multiple test adjustment}% use one of  RShowDoc("KEYWORDS")
\keyword{ effective indepent tests }% __ONLY ONE__ keyword per line
