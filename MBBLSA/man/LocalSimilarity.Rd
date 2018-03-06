\name{LocalSimilarity}
\alias{LocalSimilarity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate Local Similarity Score of Two Time Series
}
\description{
  This function computes the local similarity score for two sequences.
}
\usage{
  LocalSimilarity(x, y, maxDelay=3, rankScale = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
    two random variables with equal length.
  }
  \item{maxDelay}{
    maximum time shift allowed in computing LS score.
  }
  \item{rankScale}{
    If TRUE, perform rankNormalization first; False, otherwise.
  }
}
\value{
  \item{scoreMax}{
    Maximum local score between x and y.
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

  # Rank normalize data
  x <- rankNormalization(data[,1])
  y <- rankNormalization(data[,2])

  # compute the local similarity score
  LocalSimilarity(x, y, maxDelay = 2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

