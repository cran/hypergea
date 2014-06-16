\name{getExpectedValues}
\alias{getExpectedValues}
\title{
Expected values of multidimensional contingency tables
}
\description{
Calculates expected values of a multidimensional contingency table
}
\usage{
getExpectedValues(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
table (matrix or array)
}
}
\value{
an array containing the expected values
}
\author{
Markus Boenn
}
\examples{
x <- array(c(1:8), dim=c(2,2,2))
getExpectedValues(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ contingency table }
\keyword{ marginal }
\keyword{ expected value }
