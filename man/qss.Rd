\name{qss}
\alias{qss}
\alias{qss1}
\alias{qss2}
\alias{triogram.fidelity}
\alias{triogram.penalty}
\title{Additive Nonparametric Terms for rqss Fitting}
\description{
In the formula specification of \code{rqss} nonparametric terms 
are specified with \code{qss}.  Both univariate and bivariate
specifications are possible, and qualitative constraints may also be specified
for the qss terms.
}
\usage{
qss(x, constraint = "N", lambda = 1, w = rep(1, length(x)))
}
\arguments{
  \item{x}{The covariate determining the nonparametric component, if
	x is a matrix with two columns then  the qss
	function will construct a penalized triogram term.}
  \item{lambda}{The smoothing parameter governing the tradeoff between
	fidelity and the penalty component for this term.  In future
	versions there should be an automatic mechanism for default
	choice of the lambdas.  For now, this is the responsibility
	of the user.}
  \item{constraint}{Optional specification of qualitative constraints 
	on the fitted univariate qss functions, take the values: "N","I","D","U","C"
	"UI","UD","CI","CD" for none, increasing, decreasing, convex,
	concave, convex and increasing, etc.  And for bivariate qss
	components can take the values "N","U","C" for none, convex, and concave.} 
  \item{w}{ weights not yet unimplemented }
}
\details{
The various pieces returned are stored in sparse matrix.csr form.
See \code{rqss} for details on how they are assembled.  To preserve the
sparsity of the design matrix the first column of each qss term is dropped.
This differs from the usual convention that would have forced qss terms
to have mean zero.  This convention has implications for prediction that need to be
recognized.  The penalty components for qss terms are based on total
variation penalization of the first derivative (and gradient, for bivariate x)
as described in the references appearing in the help for \code{rqss}.
}
\value{
  \item{F}{Fidelity component of the design matrix}
  \item{A}{Penalty component of the design matrix}
  \item{R}{Constraint component of the design matrix}
  \item{r}{Constraint component of the rhs}
}
\author{Roger Koenker}

\seealso{ \code{\link{rqss}}}
\keyword{smooth}
\keyword{robust}
