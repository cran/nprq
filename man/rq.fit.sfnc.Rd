\name{rq.fit.sfnc}
\alias{rq.fit.sfnc}
\title{Sparse Constrained Regression Quantile Fitting}
\description{
Fit constrained regression quantiles using a sparse implementation of the Frisch-Newton
Interior-point algorithm.
}
\usage{
rq.fit.sfnc(x, y, R, r, tau, rhs = (1-tau)*c(t(a) %*% rep(1,length(y))), 
   nsubmax, tmpmax, nnzlmax, cachsz = 64, small = 1e-08, maxiter = 100, warn.mesg = T, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{structure of the design matrix X stored in csr format}
  \item{y}{response vector}
  \item{R}{constraint matrix stored in csr format}
  \item{r}{right-hand-side of the constraint}
  \item{tau}{desired quantile}
  \item{rhs}{the right-hand-side of the dual problem; regular users shouldn't need to specify this}
  \item{nsubmax}{upper bound of the dimension of lindx in Cholesky factorization;computed automatically inside the routine if missing}
  \item{tmpmax}{upper bound of the working array in Cholesky factorization;computed automatically inside the routine if missing}
  \item{nnzlmax}{upper bound of the non-zero entries in the Cholesky factor L;computed automatically inside the routine if missing}
  \item{cachsz}{size of the cache on the machine; default to 64}
  \item{small}{convergence tolerance for the interior algorithm}
  \item{maxiter}{upper limit for the number of iterations}
  \item{warn.mesg}{flag for printing of warning message; default to TRUE}
  \item{\dots}{additional arguments}
}
\details{
This is a sparse implementation of the Frisch-Newton algorithm for 
constrained quantile regression described in Koenker and Portnoy (1996). 
The sparse matrix
linear algebra is implemented through the functions available in the R packages
\code{SparseM}.
}
\value{
  \item{coef }{Regression quantile coefficients}
  \item{ierr }{Error code for the internal Fortran routine \code{srqfn}: 
	\describe{
		\item{1:}{ insufficient work space in call to \code{extract}}
	        \item{3:}{ insufficient storage in iwork when calling ordmmd}
      		\item{4:}{ insufficient storage in iwork when calling sfinit}
		\item{5:}{ nnzl > nnzlmax when calling sfinit}
		\item{6:}{ nsub > nsubmax when calling sfinit}
		\item{7:}{ insufficient work space in iwork when calling symfct}
		\item{8:}{ inconsistancy in input when calling symfct}
		\item{9:}{ tmpsiz > tmpmax when calling symfct; increase tmpmax}
		\item{10:}{ nonpositive diagonal encountered when calling blkfct}
		\item{11:}{ insufficient work storage in tmpvec when calling blkfct}
		\item{12:}{ insufficient work storage in iwork when calling blkfct}
		\item{13:}{ nnzd > nnzdmax in e,je when calling amub}
		\item{14:}{ nnzd > nnzdmax in g,jg when calling amub}
		\item{15:}{ nnzd > nnzdmax in h,jh when calling aplb}
			}
		}
  \item{it}{Iteration count}
  \item{time}{Amount of time used in the computation}

}
\references{Koenker, R and Ng, P. (2002).  SparseM:  A Sparse Matrix Package for \R, \cr
\url{http://www.econ.uiuc.edu/~roger/research}
}
\author{Pin Ng}
\seealso{\code{rq.fit.sfn} for the unconstrained version,
\code{SparseM} for a sparse matrix package for \R
}
\examples{
n=200
p=50
X=rnorm(n*p)
X[abs(X)<2.0]=0
X=cbind(1,matrix(X,n,p))
y=.5*apply(X,1,sum)+rnorm(n)
R=rbind(diag(p+1),-diag(p+1))
r=c(rep(0,p+1),-rep(1,p+1))
rq.o=rq.fit.sfnc(as.matrix.csr(X),y,as.matrix.csr(R),r,tmpmax=floor(1e5+exp(-12.1)*(as.matrix.csr(X)@ia[p+1]-1)^2.35))
}

\keyword{regression}
