\name{MovingBlockBootstrap}
\alias{MovingBlockBootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate statistical significance of local similarity score of two time series by Moving Block Bootstrap
}
\description{
  This function computes the statistical significance of local similarity for two sequences by Moving Block Bootstrap.
}
\usage{
  MovingBlockBootstrap(x, y, numPermu=1000, maxDelay=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
    two random variables with equal length.
  }
  \item{numPermu}{
    the number of permutation test.
  }
  \item{maxDelay}{
    maximum time shift allowed in computing LS score.
  }
}
\value{
  \item{scoreMax}{
    Maximum local score between x and y.
  }
  \item{p-value}{
    statistical significance of local score between x and y.
  }
  \item{startX}{
    start position in the most similar subintervel of x.
  }
  \item{startY}{
    start position in the most similar subintervel of y.
  }
  \item{delay}{
    time delay between x and y when the maximum local similarity score is obtained. If delay is zero, no one leads anothe; if delay is positive, y leads x; otherwise, x leads y.
  }
  \item{length}{
    length of the most similar subinterval of x and y.
  }
  \item{PosOrNeg}{
    positive of negtive relationship of x and y when the max local similarity score is obtained. If PosOrNeg is 1, positive; otherwise negative.
  }
}
\examples{
  
  # Generate two random variables
  library(MASS)
  data <- mvrnorm(100, c(0, 0), Sigma = matrix(c(1, 0.6, 0.6, 1), 2, 2))
  
  x <- data[,1]
  y <- data[,2]
  
  # compute the significance of local similarity score
  MovingBlockBootstrap(x, y, maxDelay = 2)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

