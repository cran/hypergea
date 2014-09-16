\name{NEWS}
\title{News for Package 'hypergea'}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}


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
