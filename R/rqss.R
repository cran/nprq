"[.terms" <- function (termobj, i) 
{
    resp <- if (attr(termobj, "response")) 
        termobj[[2]]
    else NULL
    newformula <- attr(termobj, "term.labels")[i]
    if (length(newformula) == 0) 
        newformula <- 1
    newformula <- reformulate(newformula, resp)
    environment(newformula) <- environment(termobj)
    terms(newformula, specials = names(attr(termobj, "specials")))
}
"rqss" <-
function (formula, tau = 0.5, data = parent.frame(), weights, 
    na.action, method = "fn", contrasts = NULL, ...) 
{
    call <- match.call()
    m <- match.call(expand = FALSE)
    temp <- c("", "formula", "data", "weights", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    m[[1]] <- as.name("model.frame")
    special <- "qss"
    Terms <- if (missing(data)) 
        terms(formula, special)
    else terms(formula, special, data = data)
    qssterms <- attr(Terms, "specials")$qss
    dropx <- NULL
    if (length(qssterms)) {
        tmpc <- untangle.specials(Terms, "qss")
        ord <- attr(Terms, "order")[tmpc$terms]
        if (any(ord > 1)) 
            stop("qss can not be used in an interaction")
        dropx <- tmpc$terms
        if (length(dropx)) 
            Terms <- Terms[-dropx]
        attr(Terms, "specials") <- tmpc$vars
    }
    m$formula <- Terms
    m <- eval(m, parent.frame())
    weights <- model.extract(m, weights)
    process <- (tau < 0 || tau > 1)
    Y <- model.extract(m, response)
    X <- model.matrix(Terms, m, contrasts)
    p <- ncol(X)
    if (length(qssterms) > 0) {
        F <- as.matrix.csr(X)
        qss <- lapply(tmpc$vars, function(u) eval(parse(text = u), 
            data))
        mqss <- length(qss)
        ncA <- rep(0, mqss + 1)
        nrA <- rep(0, mqss + 1)
        nrR <- rep(0, mqss + 1)
        for (i in 1:mqss) {
            F <- cbind(F, qss[[i]]$F)
            ncA[i + 1] <- ncol(qss[[i]]$A)
            nrA[i + 1] <- nrow(qss[[i]]$A)
            nrR[i + 1] <- ifelse(is.null(nrow(qss[[i]]$R)), 0, 
                nrow(qss[[i]]$R))
        }
        A <- as.matrix.csr(0, sum(nrA), sum(ncA))
        if (sum(nrR) > 0) {
            R <- as.matrix.csr(0, sum(nrR), sum(ncA))
            nrR <- cumsum(nrR)
        }
        ncA <- cumsum(ncA)
        nrA <- cumsum(nrA)
        for (i in 1:mqss) {
            A[(1 + nrA[i]):nrA[i + 1], (1 + ncA[i]):ncA[i + 1]] <- qss[[i]]$A
            if (nrR[i] < nrR[i + 1]) 
                R[(1 + nrR[i]):nrR[i + 1], (1 + ncA[i]):ncA[i + 
                  1]] <- qss[[i]]$R
        }
        A <- cbind(as.matrix.csr(0, nrA[mqss + 1], p), 
            A)
        if (nrR[mqss + 1] > 0) {
            R <- cbind(as.matrix.csr(0, nrR[mqss + 1], p), 
                R)
            r <- rep(0, nrR[mqss + 1])
        }
        else {
            R <- NULL
            r <- NULL
        }
        X <- rbind(F, A)
        Y <- c(Y, rep(0, nrow(A)))
	rhs <- t(rbind((1 - tau) * F, 0.5 * A))%*%rep(1,nrow(X))
	XpX <- t(X)%*%X
	nnzdmax <- XpX@ia[length(XpX@ia)]-1
	nsubmax <- max(nnzdmax,floor(1e3 + exp(-1.6) * nnzdmax^1.2))
	nnzlmax <- floor(2e5 - 2.8*nnzdmax + 7e-4 * nnzdmax^2)
	tmpmax <- floor(1e5+exp(-12.1)*nnzdmax^2.35)
        fit <- if (length(r) > 0) 
	      rqss.fit(X, Y, tau = tau, rhs = rhs, method = "sfnc", R = R,
                r = r, nsubmax = nsubmax, nnzlmax = nnzlmax, tmpmax = tmpmax)
        else 
	      rqss.fit(X, Y, tau = tau, rhs = rhs, method = "sfn",
               nsubmax = nsubmax, nnzlmax = nnzlmax, tmpmax = tmpmax)
        for(i in 1:mqss){
		ML <- p+1+ncA[i]
                MU <- p+ncA[i+1]
                qss[[i]] <- cbind(qss[[i]]$x$x,qss[[i]]$x$y,
                        c(0,fit$coef[ML:MU]))
		}
	fit$qss <- qss
    }
    else {
        fit <- if (length(weights)) 
            rq.wfit(X, Y, tau = tau, weights, method, ...)
        else rq.fit(X, Y, tau = tau, method, ...)
    }
    fit$terms <- Terms
    fit$call <- call
    fit$tau <- tau
    attr(fit, "na.message") <- attr(m, "na.message")
    class(fit) <-  "rqss"
    fit
}
"rqss.fit" <-
function (x, y, tau = 0.5, method = "sfn", ...) 
{
    fit <- switch(method, 
	sfn = rq.fit.sfn(x, y, tau = tau, ...), 
        sfnc = rq.fit.sfnc(x, y, tau = tau, ...), {
            what <- paste("rq.fit.", method, sep = "")
            if (exists(what, mode = "function")) 
                (get(what, mode = "function"))(x, y, ...)
            else stop(paste("unimplemented method:", method))
        })
    fit$contrasts <- attr(x, "contrasts")
    fit
}
"untangle.specials" <-
function (tt, special, order = 1) 
{
    spc <- attr(tt, "specials")[[special]]
    if (length(spc) == 0) 
        return(list(vars = character(0), terms = numeric(0)))
    facs <- attr(tt, "factor")
    fname <- dimnames(facs)
    ff <- apply(facs[spc, , drop = FALSE], 2, sum)
    list(vars = (fname[[1]])[spc], terms = seq(ff)[ff & match(attr(tt, 
        "order"), order, nomatch = 0)])
}
"qss" <-
function(x, constraint = "N", lambda = 1, w=rep(1,length(x))){
   if(is.list(x)) 
        if(all(!is.na(match(c("x", "y"), names(x))))) 
		qss <- qss2(x$x,x$y,constraint=constraint, lambda=lambda, w=w)
        else stop("called with list lacking xy elements")
   else if(is.matrix(x))
	if(ncol(x)==2)
		qss <- qss2(x[,1],x[,2],constraint=constraint, lambda=lambda, w=w)
        else stop("matrix argument to qss must have two columns")
   else if(is.vector(x))
	qss <- qss1(x,constraint=constraint, lambda=lambda, w=w)
   else
        stop("x argument to qss invalid")
   qss
}
"qss2" <-
function(x, y, constraint = "N", lambda = 1, w=rep(1,length(x))){
#
# Sparse Additive Quantile Smoothing Spline Models - Bivariate (Triogram) Module
#
# This function returns a structure intended to make model.matrix for a bivariate
# nonparametric component of a model formula specified by a call to rqss().  A sparse form
# of the Frisch Newton algorithm is eventually called to compute the estimator.
# An optional convexity/concavity constraint can be specified.  If
# the formula consists of a single qss component then the estimator solves the
# following variational problem:
#
#       min sum rho_tau (y_i - g(x_i)) + lambda V(grad(g))
#
# where V(f) denotes the total variation of the function f.  The solution is a piecewise
# linear function on the Delaunay triangulation formed by the observed (x_i,y_i) points.  
# Additive models can consist
# of several components of this form plus partial linear and univariate qss components.
# To resolve the identifiability problem we delete the first column of the qss design
# components.  On return F contains the fidelity portion of the design, A the penalty
# contribution of the design, R the constraint portion, and r the rhs of the constraints.
#
# Constraints are specified by the constraint argument: 
#
#	N	none
#	U	convex
#	C	concave
#
# Author:  Roger Koenker   April 2, 2003
#
# For a prototype see triogram in ~roger/projects/tv/cobar/.RData on ysidro.
#
# Warning:   Under development...todo:
#
#       o  weights
#       o  dummy x's
#       o  tau's
#       o  lambda's
#       o  ...
#
   require(SparseM)
   require(tripack)
#
   n <- length(x)
   if(n != length(y))
          stop("xy lengths do not match")
   f <- triogram.fidelity(x,y)
   F <- f$F
   A <- triogram.penalty(f$x, f$y)

switch(constraint,
        V = {   R <- A;
                r <- rep(0,nrow(R))
                },
        C = {   R <- -A;
                r <- rep(0,nrow(R))
                },
        N = { R=NULL; r=NULL}
        )
   list(x = list(x=f$x, y=f$y), F=F[,-1], A=lambda*A[,-1], R=R[,-1], r=r)
}

"qss1" <-
function(x, constraint = "N", lambda = 1, w=rep(1,length(x))){
#
# Sparse Additive Quantile Smoothing Spline Models - Univariate Module
#
# This function returns a structure intended to make model.matrix for a univariate
# nonparametric component of a model formula specified by a call to rq().  A sparse form
# of the Frisch Newton algorithm is eventually called to compute the estimator.
# Optional monotonicity and/or convexity/concavity constraints can be specified.  If
# the formula consists of a single qss component then the estimator solves the
# following variational problem:
#
#       min sum rho_tau (y_i - g(x_i)) + lambda V(g')
#
# where V(f) denotes the total variation of the function f.  The solution is a piecewise
# linear function with "knots" at the observed x_i points.  Additive models can consist
# of several components of this form plus partial linear and triogram components.
# To resolve the identifiability problem we delete the first column of the qss design
# components.  On return F contains the fidelity portion of the design, A the penalty
# contribution of the design, R the constraint portion, and r the rhs of the constraints.
#
# Constraints are specified by the constraint argument: 
#
#	N	none
#	I	monotone increasing
#	D	monotone decreasing
#	V	convex
#	C	concave
#	CI	concave and monotone increasing
#	...	etc
#
# Author:  Roger Koenker   February 27, 2003
#
# Warning:   Under development...todo:
#
#       o  weights
#       o  dummy x's
#       o  tau's
#       o  lambda's
#       o  ...
#
# See sqss for a prototype version of this function.
#
   require(SparseM)
#
   xun <- unique(x[order(x)])
   h <- diff(xun)
   nh <- length(h)
   nx <- length(x)
   p <- nh + 1
   B <- new("matrix.csr", ra = c(rbind(-1/h,1/h)), ja = as.integer(c(rbind(1:nh,2:(nh+1)))), 
		ia = as.integer(2*(1:(nh+1))-1), dimension = as.integer(c(nh,nh+1)))
   makeD <- function(p){ new("matrix.csr", ra = c(rbind(rep(-1,(p-1)),rep(1,(p-1)))), 
		ja = as.integer(c(rbind(1:(p-1),2:p))), ia = as.integer(2*(1:p)-1), 
		dimension = as.integer(c(p-1,p)))
		}
   D <- makeD(nh)
   A <- D %*% B
   if(length(xun) == length(x)){
       F <- new("matrix.csr", ra = rep(1,nx), ja = as.integer(rank(x)),
                ia = 1:(nx+1), dimension = as.integer(c( nx, nx)))
        }
   else {
       F <- new("matrix.csr", ra = rep(1,nx), ja = as.integer(factor(x)),
                ia = 1:(nx+1), dimension = as.integer(c( nx, length(xun))))
       }
   switch(constraint,
        V = {   R <- A;
                r <- rep(0,nrow(R))
                },
        C = {   R <- -A;
                r <- rep(0,nrow(R))
                },
        I = {   R <- makeD(p)
                r <- rep(0,p-1)
                },
        D = {   R <- -makeD(p)
                r <- rep(0,p-1)
                },
        VI = {  R <- makeD(p)
		R <- rbind(R,A)
                r <- rep(0,nrow(R))
		},
        VD = {  R <- -makeD(p)
		R <- rbind(R,A)
                r <- rep(0,nrow(R))
		},
        CI = {  R <- makeD(p)
		R <- rbind(R,-A)
                r <- rep(0,nrow(R))
		},
        CD = {  R <- -makeD(p)
		R <- rbind(R,-A)
                r <- rep(0,nrow(R))
		},
	N = { R=NULL; r=NULL}
	)
   list(x = list(x=xun), F=F[,-1], A=lambda*A[,-1], R=R[,-1], r=r)
}
"plot.qss1" <-
function(x, ...)
{
#plot(x[,1],x[,2],type="n",xlab="x",ylab="fit")
lines(x[,1],x[,2], ...)
}
"plot.qss2" <-
function (x, ...)
{
require(akima)
require(tripack)

        y <- x[,2]
        z <- x[,3]
        x <- x[,1]
	par(mfrow=c(1,2))
        plot(x,y,type="n",axes=TRUE,xlab="x",ylab="y")
        contour(interp(x,y,z),add=TRUE,frame.plot=TRUE, axes=TRUE, ...)
        tri <- tri.mesh(x,y)
        convex.hull(tri,plot.it=TRUE,add=TRUE)
	persp(interp(x,y,z,),theta=-40,phi=20,xlab="x",ylab="y",zlab="z", ...)
}
"plot.rqss" <-
function(x, ...) 
{
m <- length(x$qss)
for( i in 1:m){
        if(ncol(x$qss[[i]])==3){
		x$qss[[i]][,3] <- x$coef[1] + x$qss[[i]][,3]
            	plot.qss2(x$qss[[i]], ...)
	   }
        else if(ncol(x$qss[[i]])==2){
		x$qss[[i]][,2] <- x$coef[1] + x$qss[[i]][,2]
            	plot.qss1(x$qss[[i]], ...)
	   }
        else
                stop("invalid fitted qss object")
        }
}

"triogram.fidelity" <- function (x,y,ndum=0) 
{
#Make fidelity block of the triogram design in sparse matrix.csr form
#The rather esoteric match call identifies and handles duplicated xy points
n <- length(x)
A <- as.data.frame(cbind(x,y))
dupA <- duplicated(A)
if(any(dupA)){
  x <- x[!dupA]
  y <- y[!dupA]
  J <- match(do.call("paste",c(A,"\r")),do.call("paste",c(A[!dupA,],"\r")))
  z <- new("matrix.csr",ra=rep(1,n), ja=J, ia=1:(n+1),dim=as.integer(c(n,max(J))))
  }
else
  z <- as.matrix.csr(diag(n))
#Augment with dummy vertices, if any...
if(ndum > 0){
        u <- runif(ndum); v <- runif(ndum)
        xd <- min(x) + u * (max(x)-min(x))
        yd <- min(y) + v * (max(y)-min(y))
        T <- tri.mesh(x,y)
        s <- in.convex.hull(T,xd,yd)
        x <- c(x,xd[s])
        y <- c(y,yd[s])
        ndum <- sum(s)
        zdum <- matrix.csr(0,n,ndum)
        z <- cbind(z,zdum)
        }
list(x=x,y=y,F=z)
}
"triogram.penalty" <- function (x, y, eps = .Machine$double.eps) 
{
    n <- length(x)
    tri <- tri.mesh(x, y)
    bnd <- on.convex.hull(tri,x,y)
    q <- length(tri$tlist)
    m <- 13 * n
    z <- .Fortran("penalty", as.integer(n), as.integer(m), as.integer(q), 
        as.double(x), as.double(y), as.integer(bnd),as.integer(tri$tlist), 
        as.integer(tri$tlptr), as.integer(tri$tlend), rax = double(m), jax = integer(m), 
        ned = integer(1), as.double(eps), ierr = integer(1))[c("rax", 
        "jax", "iax", "ned", "ierr")]
    if (z$ierr == 1) 
        stop("collinearity in ggap")
    nnz <- 4 * z$ned
    ra <- z$rax[1:nnz]
    ja <- z$jax[1:nnz] 
    ia <- as.integer(1 + 4 * (0:z$ned))
    dim <- as.integer(c(z$ned, n))
    new("matrix.csr",ra=ra,ja=ja,ia=ia,dimension=dim)
}


