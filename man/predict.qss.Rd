\name{predict.qss}
\alias{predict.qss1}
\alias{predict.qss2}
\title{ Predict based on nonparametric quantile regression smoothing
spline component }
\description{
Additive models for nonparametric quantile regression using total
variation penalty methods can be fit with the \code{\link{rqss}}
function.  Univarariate and bivariate components can be predicted
using these functions.
}
\usage{
predict.qss1(object, newdata, ...)
predict.qss2(object, newdata, ...)
}
\arguments{
  \item{object}{ is a fitted object produced by \code{\link{rqss}} }
  \item{newdata}{ a data frame describing the observations at which
	prediction is to be made   }
  \item{\dots}{ optional arguments }
}
\details{
For both univariate and bivariate prediction linear interpolation is
done.  In the bivariate case, this involves computing barycentric
coordinates of the new points relative to their enclosing triangles.
}
\value{
A list consisting of x and y components in the case of qss1,
and a list consisting of x, y, and z  components in the case of qss2.
In the former case the y component constitutes the predictions at x,
and in the latter, z is the vector of the predictions at the points
(x,y).
}
\author{ R. Koenker }
\seealso{  \code{\link{rqss}}}
\examples{ }
\keyword{regression}
\keyword{smooth}
\keyword{robust}
