#-------------------------------------------------------------------------------
".First.lib" <-
function(lib, pkg) {
   require(SparseM)
   library.dynam("nprq", pkg, lib)
   print("nprq loaded")
}
#-------------------------------------------------------------------------------
"rq.fit.sfn" <- function(a,y,tau=.5, rhs = (1-tau)*c(t(a) %*% rep(1,length(y))),
	nsubmax,tmpmax,nnzlmax,cachsz=64,small=1.e-6,maxiter=100,warn.mesg=TRUE,...){
#################################################################
# Interface for a sparse implementation of LMS interior point method
#	a -- structure of the design matrix A' stored in csr format
#	y -- pseudo response vector
#	m -- column dimension of A' in full mode
#       rhs -- the rhs of the dual problem -- specify at your own risk
#	tau -- the desired regression quantile
#	nnzemax -- maximum number of non-zero elements in AA'
#	tmpmax -- upper bound for the dimension of tmpvec
#	cachsz -- size of the cache memory; it's machine dependent
#	small -- the tolerance used to check for convergence
#################################################################
	y <- -y
	n <- length(y)
	m <- a@dimension[2]
	if(n != a@dimension[1]) 
		stop("Dimensions of design matrix and the response vector are not compatible")
	u <- rep(1,length=n)
	x <- rep((1-tau),length=n)
	nnzdmax <- nnza <- a@ia[n+1]-1
	iwmax <- 7*m+3
	ao <- t(a)
	e <- ao%*%a
	nnzemax <- e@ia[m+1]-1
	if(missing(nsubmax)) nsubmax <- nnzemax
	if(missing(nnzlmax)) nnzlmax <- 4*nnzdmax
	if(missing(tmpmax)) tmpmax <- 6*m
	level <- 8
	wwm <- vector("numeric",3*m)
	s <- u - x
	b1 <- solve(e,ao%*%y,tmpmax=tmpmax,nnzlmax=nnzlmax,nsubmax=nsubmax)
	r <- y - a%*%b1
	z <- ifelse(abs(r)<small,(r*(r>0)+small),r*(r>0))
	w <- z - r
	wwn <- matrix(0,n,14)
	wwn[,1] <- r
	wwn[,2] <- z
	wwn[,3] <- w
	srqfnb.o <- .Fortran("srqfn",
		n = as.integer(n),
		m = as.integer(m),
		nnza = as.integer(nnza),
		a = as.double(a@ra),
		ja = as.integer(a@ja),
		ia = as.integer(a@ia),
		ao = as.double(ao@ra),
		jao = as.integer(ao@ja),
		iao = as.integer(ao@ia),
		nnzdmax = as.integer(nnzdmax),
		d = double(nnzdmax),
		jd = integer(nnzdmax),
		id = integer(m+1),
		dsub = double(nnzemax+1),
		jdsub = integer(nnzemax+1),
		nnzemax = as.integer(nnzemax),
		e = as.double(e@ra),
		je = as.integer(e@ja),
		ie = as.integer(e@ia),
		nsubmax = as.integer(nsubmax),
		lindx = integer(nsubmax),	
		xlindx = integer(m+1),
		nnzlmax = as.integer(nnzlmax),
		lnz = double(nnzlmax),
		xlnz = integer(m+1),
		iw = integer(m*5),
		iwmax = as.integer(iwmax),
		iwork = integer(iwmax),
		xsuper = integer(m+1),
		tmpmax = as.integer(tmpmax),
		tmpvec = double(tmpmax),
		wwm = as.double(wwm),
		wwn = as.double(wwn),
		cachsz = as.integer(cachsz),
		level = as.integer(level),
		x = as.double(x),
		s = as.double(s),
		u = as.double(u),
		c = as.double(y),
		sol = as.double(b1),
		rhs = as.double(rhs),
		small = as.double(small),
		ierr = integer(1),
		maxiter = as.integer(maxiter),
		time = double(7),
		PACKAGE = "nprq")
	if(!all(srqfnb.o$ierr==0) && warn.mesg) warning("Check ierr!\n")
	list(coef=-srqfnb.o$sol,ierr=srqfnb.o$ierr,it=srqfnb.o$maxiter,time=sum(srqfnb.o$time))
}
#-------------------------------------------------------------------------------
"rq.fit.sfnc" <- function(x,y,R,r,tau=0.5, rhs = (1-tau)*c(t(x) %*% rep(1,length(y))),
	nsubmax,tmpmax,nnzlmax,cachsz=64,small=1.e-8,maxiter=100,warn.mesg=TRUE,...){
#################################################################
# Interface for a sparse implementation of LMS interior point method
#	x -- structure of A1' stored in csr format
#	y -- pseudo response vector
#	m -- column dimension of A1' in full mode
#       R -- structure of the constraint matrix A2' stored in csr format
#	r -- rhs of the inequality constraints
#	tau -- the desired regression quantile
#       rhs -- the rhs of the dual problem -- specify at your own risk
#	nnzhmax -- upper bound of the non-zero elements in AA'
#       cachsz -- size of the cache memory; it's machine dependent
#       small -- the tolerance used to check for convergence
#################################################################
	y <- -y
	r <- -r
	n1 <- length(y)
	m <- x@dimension[2]
	if(n1 != x@dimension[1]) stop("The design matrix A1' and response vector y are not compatible")
	n2 <- length(r)
	if(n2 != R@dimension[1]) stop("The constraint matrix A2' and the constraint right-hand-side are not compatible")
	maxn1n2 <- max(n1,n2)
	u <- rep(1,length=n1)
	x1 <- rep(1-tau,length=n1)
	x2 <- rep(1,length=n2)
	wwm <- vector("numeric",6*m)
	wwm[1:m] <- rhs
	nnzx <- x@ia[x@dimension[1]+1]-1
	nnzR <- R@ia[R@dimension[1]+1]-1
	nnzdmax <- max(nnzx,nnzR)
	iwmax <- 7*m+3
	ao1 <- t(x)
	ao2 <- t(R)
	e <- ao1%*%x
	g <- ao2%*%R
	h <- e + g
	nnzemax <- e@ia[e@dimension[1]+1]-1
	nnzgmax <- g@ia[g@dimension[1]+1]-1
	nnzhmax <- h@ia[h@dimension[1]+1]-1
	if(missing(nnzlmax)) nnzlmax <- 4*nnzdmax 
	if(missing(nsubmax)) nsubmax <- nnzhmax
	if(missing(tmpmax)) tmpmax <- 6*m
	s <- u - x1
	chol.o <- chol(ao1%*%x,tmpmax=tmpmax,nsubmax=nsubmax,nnzlmax=nnzlmax)
	b <- backsolve(chol.o,ao1%*%y)
	r1 <- y - x%*%b
	z1 <- ifelse(abs(r1)<small,(r1*(r1>0)+small),r1*(r1>0))
	w <- z1 - r1
	z2 <- rep(1,n2)
	wwn1 <- matrix(0,n1,10)
	wwn1[,1] <- z1
	wwn1[,2] <- w
	wwn2 <- matrix(0,n2,7)
	wwn2[,2] <- z2
	level <- 8
	srqfnc.o <- .Fortran("srqfnc",
		n1 = as.integer(n1),
		m = as.integer(m),
		nnzx = as.integer(nnzx),
		x = as.double(x@ra),
		jx = as.integer(x@ja),
		ix = as.integer(x@ia),
		ao1 = as.double(ao1@ra),
		jao1 = as.integer(ao1@ja),
		iao1 = as.integer(ao1@ia),
		n2 = as.integer(n2),
		nnzR = as.integer(nnzR),
		R = as.double(R@ra),
		jR = as.integer(R@ja),
		iR = as.integer(R@ia),
		ao2 = as.double(ao2@ra),
		jao2 = as.integer(ao2@ja),
		iao2 = as.integer(ao2@ia),
		nnzdmax = as.integer(nnzdmax),
		d = double(nnzdmax),
		jd = integer(nnzdmax),
		id = integer(m+1),
		dsub = double(nnzhmax+1),
		jdsub = integer(nnzhmax+1),
		nnzemax = as.integer(nnzemax),
		e = as.double(e@ra),
		je = as.integer(e@ja),
		ie = as.integer(e@ia),
		nnzgmax = as.integer(nnzgmax),
		g = double(nnzgmax),
		jg = integer(nnzgmax),
		ig = integer(m+1),
		nnzhmax = as.integer(nnzhmax),
		h = double(nnzhmax),
		jh = integer(nnzhmax),
		ih = integer(m+1),
		nsubmax = as.integer(nsubmax),
		lindx = integer(nsubmax),	
		xlindx = integer(m+1),
		nnzlmax = as.integer(nnzlmax),
		lnz = double(nnzlmax),
		xlnz = integer(m+1),
		iw = integer(m*5),
		iwmax = as.integer(iwmax),
		iwork = integer(iwmax),
		xsuper = integer(m+1),
		tmpmax = as.integer(tmpmax),
		tmpvec = double(tmpmax),
		maxn1n2 = as.integer(maxn1n2),
		ww1 = double(maxn1n2),
		wwm = as.double(wwm),
		wwn1 = as.double(wwn1),
		wwn2 = as.double(wwn2),
		cachsz = as.integer(cachsz),
		level = as.integer(level),
		x1 = as.double(x1),
		x2 = as.double(x2),
		s = as.double(s),
		u = as.double(u),
		c1 = as.double(y),
		c2 = as.double(r),
		sol = as.double(b),
		small = as.double(small),
		ierr = integer(1),
		maxiter = as.integer(maxiter),
		time = double(7),
		PACKAGE = "nprq")
	if(srqfnc.o$ierr==13 && warn.mesg) stop("Increase nnzh.factor")
	if(!all(srqfnc.o$ierr==0) && warn.mesg) warning("Check ierr!\n")
	list(coef=-srqfnc.o$sol,ierr=srqfnc.o$ierr,it=srqfnc.o$maxiter,time=sum(srqfnc.o$time))
}
