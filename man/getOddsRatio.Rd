\name{getOddsRatio}
\alias{getOddsRatio}
\title{
Odds ratio
}
\description{
This function calculates odds ratios for 2x2 and 2x2x2 contingency tables.
}
\usage{
getOddsRatio(x, za=TRUE)
}

\arguments{
  \item{x}{
a contingency table (array or matrix)
}
  \item{za}{
zero adjustment, i.e. how to behave in case of zeros (logical). If FALSE than odds ratios as calculated. If TRUE, each cell is adjusted by addition of a small number (0.5).
}

}
\details{
If \code{x} is a 2x2 table then the odds ratio is returned. If \code{x} is a 2x2x2 table then Bartlett's ratio of odds ratios is calculated.

To be able to calculate odds ratios even for huge numbers, the log of each cell count is taken.

Calculation of an odds ratio becomes difficult in presence of zeros. Depending on the position(s) of zero(s) in the table, the ratio becomes zero, NaN or Inf. This corresponds to the behaviour of the \code{\link{fisher.test}} in such a case. However, another strategy is chosen by the vcd pacakge. Here a small number (0.5) is added to each cell.

Note that this routine calculates the sample odds ratio (also: conditional Maximum Likelihood estimate). This is different to other routines like \code{\link{fisher.test}}, which calculate the conditional Maximum Likelihood Estimate (MLE).

}
\value{
the odds ratio (numeric)
}
\references{
    MARTIN A. HAMILTON (1979) CHOOSING THE PARAMETER FOR A 2 x 2 TABLE OR A 2 x 2 x 2 TABLE ANALYSIS. Am. J. Epidemiol. (1979) 109 (3): 362-375
}
\author{
Markus Boenn
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

%Cinfidence interval fuer odds ratio
%http://www.biochemia-medica.com/content/odds-ratio-calculation-usage-and-interpretation
%odds ratio bei fishertest
%http://www.utstat.toronto.edu/~brunner/oldclass/312f12/lectures/312f12FisherWithR.pdf
%pruefe ob <environment: namespace:stats> bei hypergea alles analog so steht