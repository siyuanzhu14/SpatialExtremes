
library(SpatialExtremes)
data(rainfall)

data(swissalt)##elevation data for Switzerland

jpeg("swissstations.jpeg")
par(mar = rep(0, 4), ps = 16)
image(lon.vec, lat.vec, alt.mat, col = terrain.colors(64), asp = 1,bty = "n", xlab = "n", ylab = "n", axes = FALSE)
swiss(add = TRUE, city = TRUE)
points(coord, pch = 15)
dev.off()


sites = c(3,4,5,6)
loc.form <- scale.form <- y ~ lon + lat; shape.form <- y ~ 1

## Transform rain data to to unit Fréchet first
frechet <- apply(rain, 2, gev2frech, emp = TRUE)

fit_copula_gaussian_exp <- fitcopula(frechet, coord[,1:2], "gaussian","powexp", nugget=0, loc.form, scale.form, shape.form)

fit_copula_gaussian_whitmat <- fitcopula(frechet, coord[,1:2], "gaussian","whitmat", nugget=0, loc.form, scale.form, shape.form)

fit_copula_student_whitmat <- fitcopula(frechet, coord[,1:2], "student","whitmat", nugget=0, loc.form, scale.form, shape.form)

fit_copula_student_exp <- fitcopula(frechet, coord[,1:2], "student","powexp", nugget=0, loc.form, scale.form, shape.form)



#schlather
fit_whitmat <- fitmaxstab(frechet, coord[,1:2], "whitmat", nugget = 0)

fit_cauchy <- fitmaxstab(frechet, coord[,1:2], "cauchy", nugget = 0)

#geometric gaussian whith whitmat
fit_gwhitmat <- fitmaxstab(frechet, coord[,1:2], "gwhitmat", nugget = 0)

#smith
fit_gauss <- fitmaxstab(frechet, coord[,1:2], "gauss", iso = TRUE)

fit_brown <- fitmaxstab(frechet, coord[,1:2], "brown")

fit_exp <- fitmaxstab(frechet, coord[,1:2], "powexp", nugget = 0)







#latent variable
hyper <- list()
hyper$betaMeans <- list(loc = rep(0, 3), scale = rep(0, 3), shape = 0)
hyper$betaIcov <- list(loc = diag(rep(1/10, 3)), scale = diag(rep(1/10, 3)), shape = 1/10)
#………………
hyper$sills <- list(loc = c(1, 12), scale = c(1, 1), shape = c(1, 0.04))
hyper$ranges <- list(loc = c(5, 3), scale = c(5, 3), shape = c(5, 3))
hyper$smooths <- list(loc = c(1, 1), scale = c(1, 1), shape = c(1, 1))
#………………………………………………………………………………………………………………………………………………………………
prop <- list(gev = c(3, 0.1, 0.3), ranges = c(1, 0.8, 1.2), smooths = rep(0, 3))
start <- list(sills = c(10, 10, 0.5), ranges = c(20, 10, 10), smooths = c(1, 1, 1), beta = list(loc = c(25, 0, 0), scale = c(33, 0, 0), shape = 0.001)) 
#OTLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL	
res <- latent(rain, coord[,1:2], "powexp", loc.form, scale.form, shape.form, hyper = hyper, prop = prop,start = start, n = 300000, burn.in = 5000, thin = 30); 
res




#???????????????????????????????????????
variogram(rain, coord[,1:2], col = "grey")
variogram(data = fit_copula_gaussian_exp,coord = coord[,1:2], n.bins = 200, add = TRUE)
#??????????????????????????????????????

jpeg("fig6_2.jpeg")
fmadogram(rain, coord[,1:2], which = "ext", col = "grey", main = "Husler Reiss")
fmadogram(fitted = fit_gwhitmat, which = "ext", n.bins = 200, add = TRUE)
plot(fit_copula_gaussian_exp$ext.coeff, from = 0, to = 110, col = 1, add = TRUE)
plot(fit_copula_student_exp$ext.coeff, from = 0, to = 110, col = 3, add = TRUE)
legend("bottomright", c("Gauss Copula", "Stable", "Exponential"),
       col = 1:3, lty = 1, bty = "n", inset = 0.05)
dev.off()

jpeg("fig6_3.jpeg")
fmadogram(rain, coord[,1:2], which = "ext", col = "grey", main = "Extremal t")
fmadogram(fitted = fit_gwhitmat, which = "ext", n.bins = 200, add = TRUE)
plot(fit_copula_student_whitmat$ext.coeff, from = 0, to = 110, col = 1, add = TRUE)
plot(fit_cauchy$ext.coeff, from = 0, to = 110, col = 3, add = TRUE)
legend("bottomright", c("Whittle-St Copula", "Whittle", "Cauchy"),
       col = 1:3, lty = 1, bty = "n", inset = 0.05)
dev.off()



jpeg("fig9_1.jpeg")
fmadogram(rain, coord[,1:2], which = "ext", col = "grey",main = "Schlather")
fmadogram(fitted = fit_whitmat, which = "ext", n.bins = 200, add = TRUE)## <<- you can pass a fitted model also
plot(fit_cauchy$ext.coeff, from = 0, to = 110, col = 3, add = TRUE)
legend("bottomright", c("Whittle", "Cauchy"),
       col = 2:3, lty = 1, bty = "n", inset = 0.05)
dev.off()

jpeg("fig9_2.jpeg")
fmadogram(rain, coord[,1:2], which = "ext", col = "grey", main = "Geometric")
fmadogram(fitted = fit_gwhitmat, which = "ext", n.bins = 200, add = TRUE)## <<- you can pass a fitted model also
plot(fit_cauchy$ext.coeff, from = 0, to = 110, col = 3, add = TRUE)
legend("bottomright", c("Whittle", "Cauchy"),
       col = 2:3, lty = 1, bty = "n", inset = 0.05)
dev.off()

jpeg("fig9_3.jpeg")
fmadogram(rain, coord[,1:2], which = "ext", col = "grey", main = "Brown-Resnick")
fmadogram(fitted = fit_whitmat, which = "ext", n.bins = 200, add = TRUE)## <<- you can pass a fitted model also
plot(fit_gauss$ext.coeff, from = 0, to = 110, col = 3, add = TRUE)
legend("bottomright", c("Brown-Resnick", "Smith"),
       col = 2:3, lty = 1, bty = "n", inset = 0.05)
dev.off()


