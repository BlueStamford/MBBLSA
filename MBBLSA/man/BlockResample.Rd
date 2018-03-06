\name{BlockResample}
\alias{BlockResample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Moving block bootstrap resample
}
\description{
  This function achieves the moving block bootstrap resample of sequence.
}
\usage{
  BlockResample(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    two random variables with equal length.
  }
}
\value{
This function returns the resampled sequence of x by moving block bootstrap.
}
\examples{
  
  # Generate random variables
  x <- rnorm(100)
  BlockResample(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

