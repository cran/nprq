
#Demo for an rgl Animation of Cobar Ore fitting 

#Make an initial quite rough fit of the data
library(nprq)
library(rgl)
data(CobarOre)
fit <- rqss(z ~ qss(cbind(x,y), lambda = .01, ndum=100),data =
CobarOre)
dummies <- fit$qss[[1]]$dummies
zcol <- CobarOre$z
plot(fit,rgl = TRUE) 
cat("Now orient the plot as needed: \n Resize window, \n mouse button 1 to change viewpoint, \n mouse button 2 to zoom, \n and hit return when ready") 
scan()
rgl.bg(color="8")
for(i in 1:20){
        fname <- paste("cobar",i,".png",sep="")
        lam <- 2*i/100 
        fit <- rqss(z ~ qss(cbind(x,y), lambda = lam, dummies = dummies),data = CobarOre)
        rgl.clear()
        plot(fit,rgl = TRUE, zcol = zcol) 
        rgl.snapshot(fname)
        }       
