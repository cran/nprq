\name{plot}
\alias{plot.rqss}
\alias{plot.qss1}
\alias{plot.qss2}
\title{Default Ploting Method for rqss()}
\description{
Takes a fitted \code{rqss} object produced by \code{rqss()} and plots
the component smooth functions that make up the ANOVA decomposition.
Since the components "omit the intercept" the estimated intercept is added back 
in -- this facilitates the comparison of quantile fits particularly.
}
\usage{
plot.rqss(x, ...)
}
\arguments{
  \item{x}{a fitted \code{rqss} object produced by \code{rqss()}.}
  \item{...}{
    additional arguments for the plotting algorithm}
}
\value{
    The function produces plots for the ANOVA components as a side effect.
}
\author{ Roger Koenker }
\seealso{ \code{\link{rqss}}}
\examples{
n <- 200
x <- sort(rchisq(n,4))
z <- x + rnorm(n)
y <- log(x)+ .1*(log(x))^2 + log(x)*rnorm(n)/4 + z
plot(x,y-z)
fit <- rqss(y~qss(x,constraint="N")+z)
lines(x[-1],fit$coef[1]+fit$coef[-(1:2)])
fit <- rqss(y~qss(x,constraint="I")+z)
lines(x[-1],fit$coef[1]+fit$coef[-(1:2)],col="blue")
fit <- rqss(y~qss(x,constraint="CI")+z)
lines(x[-1],fit$coef[1]+fit$coef[-(1:2)],col="red")
#Cleanup
rm(list=ls())
#A bivariate example
data(CobarOre)
attach(CobarOre)
fit <- rqss(z~qss(cbind(x,y),lambda=.08))
plot(fit)
}
\keyword{regression}
\keyword{smooth}
\keyword{iplot}
