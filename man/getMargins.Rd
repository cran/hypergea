\name{getMargins}
\alias{getMargins}
\title{
Marginal sums of multidimensional contingency tables
}
\description{
Calculates marginal sums for IxJ and IxJxK tables
}
\usage{
getMargins(x)
}
\arguments{
  \item{x}{
a table (matrix, array) without missing values
}
}
\value{
A list of length \code{length(dim(x))} containing the marginal sums of each dimension
}
\examples{
x <- matrix(c(1:4), nrow=2, ncol=2)
getMargins(x)

x <- matrix(c(1:6), nrow=2, ncol=3)
getMargins(x)

x <- array(c(1:8), dim=c(2,2,2))
getMargins(x)
}
