\name{getOddsRatio}
\alias{getOddsRatio}
\title{
Odds ratio
}
\description{
This function calculates odds ratios for \eqn{2\times2}{2x2} and \eqn{2\times2\times2}{2x2x2} contingency tables.
}
\usage{
getOddsRatio(x, za=TRUE)
}

\arguments{
  \item{x}{
a contingency table (array or matrix)
}
  \item{za}{
zero adjustment, i.e. how to behave in case of zeros (logical). If \code{FALSE} than odds ratios as calculated. If \code{TRUE}, each cell is adjusted by addition of a small number (0.5).
}

}
\details{
If \code{x} is a \eqn{2\times2}{2x2} table then the odds ratio is returned, defined as
\deqn{ OR=\frac{x_{00} x_{11} }{x_{10} x_{01}} }{ OR=(x_00 x_11 )/( x_10 x_01 ) }
If \code{x} is a \eqn{2\times2\times2}{2x2} table then Bartlett's ratio of odds ratios is calculated, defined as
\deqn{ OR=\frac{x_{000} x_{111} x_{011} x_{101} }{ x_{101} x_{001} x_{100} x_{111} } }{ OR=(x_000 x_111 x_011 x_101 )/( x_101 x_001 x_100 x_111 ) }
To be able to calculate odds ratios even for huge numbers, the log of each cell count is taken.

Calculation of an odds ratio becomes difficult in presence of zeros. Depending on the position(s) of zero(s) in the table, the ratio becomes zero, \code{NaN} or \code{Inf}. This corresponds to the behaviour of the \code{\link{fisher.test}} in such a case. However, another strategy is (Haldane) to add a small number (\eqn{0.5}) to each cell when the table contains at least one zero.

Note that this routine calculates the sample odds ratio (also: unconditional Maximum Likelihood estimate). This is different to other routines like \code{\link{fisher.test}}, which calculate the conditional Maximum Likelihood Estimate (MLE).

}
\value{
the odds ratio (numeric) if \code{x} is a \eqn{2\times2}{2x2} or \eqn{2\times2\times2}{2x2x2} table, \code{NA} otherwise
}
\references{
    MARTIN A. HAMILTON (1979) CHOOSING THE PARAMETER FOR A 2 x 2 TABLE OR A 2 x 2 x 2 TABLE ANALYSIS. Am. J. Epidemiol. (1979) 109 (3): 362-375
}

\seealso{
\code{\link{fisher.test}}, \code{vcd} package for visualisation of contingency tables
}

\examples{
CT <- matrix(c(1:4), nrow=2)
getOddsRatio(CT)

CT <- array(c(1:8), dim=c(2,2,2))
getOddsRatio(CT)
}
\keyword{ contingency table }
\keyword{ odds ratio }

