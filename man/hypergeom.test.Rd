\name{hypergeom.test}
\alias{hypergeom.test}
\title{
Hypergeometric test
}
\description{
Performs exact hypergeometric tests for \eqn{I \times J} and \eqn{2 \times 2 \times 2} contingency tables
}
\usage{
hypergeom.test(x, alternative = "two.sided",
	pval.method = "minimum.likelihood", nthreads = 2, \dots)
}
\arguments{
  \item{x}{
an integer array (or matrix) without missing values
}
  \item{alternative}{
a character string specifying the alternative hypothesis,
          must be one of \dQuote{two.sided} (default), \dQuote{greater} or
          \dQuote{less}.
}
  \item{pval.method}{
a character string specifying which approach should be used to determine the p-value in a two-sided test. Can be \dQuote{minimum.likelihood} or \dQuote{double}. Ignored if \sQuote{alternative} is \dQuote{less} or \dQuote{greater}
}
  \item{nthreads}{
the number of threads that should be used
}
  \item{...}{
additional arguments. Currently not used
}

}
\details{
There is a lot of discussion about the correct way to calculate a p-value in a two sided test. Thus, several alternatives are possible (see Agresti (1992) for a short discussion and further references).

Currently two types of two-sided p-values are implemented. The first one (default) is the minimum-likelihood approach, which is also used by \code{\link{fisher.test}}, for instance. The second one is the so-called doubling approach.

}
\value{
A list with class \sQuote{"htest"} containing the following components
\item{statistic}{the value observed if the cell of interest}
\item{p.value}{the p-value for the selected test}
\item{estimate}{the odds ratio calculated from \sQuote{x}}
\item{alternative}{a character string describing the alternative hypothesis}
\item{method}{a character string indicating what type of hypergeometric test was performed}
\item{conf.int}{the 95\% confidence interval of the odds ratio}
}

\source{
The code for 2x2x2 hypergeometric test is a C translation of the
     Fortran algorithm of Zelterman (1999).
}

\references{
	Agresti, Alan (1992) A survey of exact inference for contingency tables. Statistical Science Vol. 7, No. 1 (Feb., 1992), pp. 131-153

	Zelterman, Daniel (1999) Models for Discrete Data. Oxford University Press
}

\seealso{
\code{\link{fisher.test}} to calculate p-values for \eqn{2\times 2} tables
}
\examples{
	##########
	## example from Zelterman (1999), p. 119, Table 4.1
	## mice were exposed to a fungicide (or not), some have a tumor
	##########
	dmn <- list(Tumor=c('y', 'n'), Exposition=c('y', 'n'))
	CT <- matrix(c(4,12,5,74), nrow=2, dimnames=dmn)
	\dontrun{CT}
	hypergeom.test(CT)

	\dontrun{
	##########
	## example from Agresti (1992), Table 2
	## Fisher's tea tasting experiment
	# Chi-squared based tests are underestimating p-value
	##########
	dmn <- list('poured first'=c('m', 't'), guess=c('m', 't'))
	CT <- matrix(c(3,1,1,3), nrow=2, dimnames=dmn)
	hypergeom.test(CT, alternative="two.sided")
	hypergeom.test(CT, alternative="two.sided", pval.method="double")
	hypergeom.test(CT, alternative="less")
	hypergeom.test(CT, alternative="greater")
	}


	##########
	## example from Zelterman (1999), p. 44, Table 1.2
	## Response to three drugs
	##########
	dmn <- list(A=c('F', 'U'), B=c('F', 'U'), C=c('F', 'U'))
	CT <- array(c(6,16,2,4,2,4,6,6), dim=c(2,2,2), dimnames=dmn)
	\dontrun{CT}
	hypergeom.test(CT)

}
\keyword{ hypergeometric }
\keyword{ test }
\keyword{ OpenMP }
