if(1==2){

library(PharmacoGx)
GDSC <- readRDS("~/CXP/CL/Data/GDSC_30Nov2015.Rda")
#GDSC@sensitivity



z=dc.fit[[2]]
par(pty="s")
plot(z$data$x, z$data$y)
abline(z$fit)
Point0 = z$data$y[1]

lxT = predict(z$fit, newdata = list(x=0))+ Point0
abline(lxT, coef(z$fit), col="#a50f15")

x = z$data$x
y = (coef(z$fit)* x) + lxT
lines(x, y)


ang = atan(coef(z$fit)[["I(x - p)"]])
yt = (tan(ang)* x) + lxT
lines(x, yt)


#
}


