\name{NEWS}
\title{News for Package 'hypergea'}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}


\section{Changes in hypergea version 1.3.4 (2017-????)}{
  \itemize{
    \item Ensure that P-value lies between zero and one, which is not always true for 'double' approach applied on very sparse matrices with row-/col-sums equal to zero.
  }
}


\section{Changes in hypergea version 1.3.3 (2016-09-09)}{
  \itemize{
    \item \code{hypergeom.test} now throws warning in case of p-value identical to zero (two-sided test) and changes p-value using other p-value-approach
    \item contact e-mail address changed to \code{markus.boenn.sf@googlemail.com}
    \item minor changes in C-code \code{run2x2} and \code{runIx2}
    \item unit tests added
    \item added functions for exclusive use by unit tests
  }
}


\section{Changes in hypergea version 1.3.2 (2016-04-14)}{
  \itemize{
    \item added CITATION-file
  }
}


\section{Changes in hypergea version 1.3.1 (2016-04-09)}{
  \itemize{
    \item fixed 'bug' in source \code{hypergeom_IxJ.cpp} causing compile error on Solaris (relevant for R 3.3.0 and higher)
  }
}


\section{Changes in hypergea version 1.3.0 (2015-01-02)}{
  \itemize{
    \item added \code{hypergeom.test.list}, which is performs \code{hypergeom.test} to a list of contingency tables with same population size. Only available for 2x2 tables until full code optimization.
    \item fixed typo in output formation of test result
  }
}


\section{Changes in hypergea version 1.2.3 (2014-09-16)}{
  \itemize{
    \item in case of two-sided test, the original minimum-likelihood approach for p-value calculation is not the default anymore. Instead the R-implementation in \code{fisher.test} is used.
    \item minor change in C-routine for 2x2x2 tables
  }
}


\section{Changes in hypergea version 1.2.2 (2014-09-12)}{
  \itemize{
    \item bug-fix in routine for 2x2 tables
     \item IxJ table transposed if J>I for speed improvement
  }
}


\section{Changes in hypergea version 1.2.1 (2014-09-09)}{
  \itemize{
    \item specialized C-routines for 2x2, Ix2, IxJ contingency tables for speed improvement
    \item 2xI table transposed to Ix2 table for speed improvement
  }
}
