"predict.qss1" <-
function(object,  newdata, ...){
	x <- object[,1]
	y <- object[,2]
	if(any(diff(x)<0))stop("x coordinates in object not monotone")
	if(max(newdata) > max(x) || min(newdata) < min(x))
		stop("no extrapolation allowed in predict.qss")
	bin  <- cut(newdata,x,label=FALSE) 
	p <- length(x)
	m <- length(newdata)
	V <- cbind(bin,bin+1)
	B <- cbind(x[bin+1]-newdata, newdata - x[bin])/(x[bin+1]-x[bin]) 
	D <- new("matrix.csr",ra = c(t(B)),ja=as.integer(c(t(V))),
		ia = as.integer(c(2*(1:m)-1,2*m+1)),dim=c(m,p))
	list(x = newdata,y = D%*%y)
        }
"predict.qss2" <-
function(object, newdata,  ...){
        #predict at points in newdata using the object (from rqss)
	x <- object[,1]
	y <- object[,2]
	z <- object[,3]
	tri.area <- function(v){
	  0.5*((v[2,1]-v[1,1])*(v[3,2]-v[1,2])-(v[3,1]-v[1,1])*(v[2,2]-v[1,2]))
	  }
	barycentric <- function(v){
	  b <- rep(0,3)
	  Area <- tri.area(v[1:3,])
	  b[1] <- tri.area(v[c(4,2,3),])/Area
	  b[2] <- tri.area(v[c(1,4,3),])/Area
	  b[3] <- tri.area(v[c(1,2,4),])/Area
	  if(any(b < 0 || b > 1))stop("barycentric snafu")
	  b
	  }
        if(is.list(newdata))
                if(all(!is.na(match(c("x","y"),names(newdata))))){
                        newx <- newdata$x
                        newy <- newdata$y
                        }
                else(stop("newdata list must have x and y elements"))
        else if(is.matrix(newdata))
                if(ncol(newdata)==2){
                        newx <- newdata[,1]
                        newy <- newdata[,2]
                        }
                else(stop("newdata matrix must have 2 columns"))

        tri <- tri.mesh(x,y)
        if(!all(in.convex.hull(tri,newx,newy)))
                stop("some newdata points outside convex hull")
        p <- length(x)
        m <- length(newx)
        V <- matrix(0,m,3)
        B <- matrix(0,m,3)
        for(i in 1:m){
                V[i,] <- unlist(tri.find(tri,x[i],y[i]))
                v <- rbind(cbind(object$x[V[i,]],object$y[V[i,]]),c(x[i],y[i]))
                B[i,] <- barycentric(v)
                }
	D <- new("matrix.csr",ra = c(t(B)),ja=as.integer(c(t(V))),ia =
		as.integer(c(4*(1:m)-1,4*m+1)),dim=c(m,p))
	list(x=newx,y=newy,z=D%*%z)
        }
