"lprq" <-
function(x, y, h, tau = .5, m = 50)
{
	require(quantreg)
        xx <- seq(min(x),max(x),length=m)
        fv <- xx
        der <- xx
        for(i in 1:length(xx)) {
                z <- x - xx[i]
                wx <- dnorm(z/h)
                r <- rq(y~z, weights=wx, tau=tau, ci=FALSE)
                fv[i] <- r$coef[1.]
                der[i] <- r$coef[2.]
        }
        list(xx = xx, fv = fv, der = der)
}

