stop("CODE HAS BEEN ABANDONED, THE RVM WAS REMOVED FROM FURTHER STUDY, ",
     "AND ProbSVMexp.R IS CANONICAL TO THE PAPER.")

library(kernlab)
library(mgcv)
seed <- 1

NSE <- function(obs,sim) {
   1 - (sum((obs-sim)^2)/sum((obs-mean(obs))^2))
}

par(las=2)

xymod <- function(x,y, datum=70) {
  i <- 20*sin(+1.75*pi*(x-y)/100)
  j <- 20*cos(-1.50*pi*   y /100)
  z <- i + j
  err <- rnorm(length(x), mean=0, sd=2)
  drift <- 0.001*x^2 - 0.1*y
  z <- z + drift + err + datum
  return(z)
}

grid <- 1:200
M <- matrix(nrow=length(grid), ncol=length(grid))
for(i in 1:length(grid)) {
  for(j in 1:length(grid)) {
    M[i,j] <-  xymod(grid[i], grid[j])
  }
}

set.seed(62)
nsim <- 800
X <- runif(nsim, min=0, 200)
Y <- runif(nsim, min=0, 200)
Z <- xymod(X,Y)

set.seed(seed)

#pdf("../draftfigures/RVMNOTUSEDsurface.pdf", useDingbats=FALSE)
plot(X/200, Y/200, type="n",
     xlab="FRACTION EASTING RANGE", ylab="FRACTION NORTHING RANGE")
contour(M, lwd=1.5,labcex = 0.7, nlevels=8, add=TRUE, axes=FALSE)
points(X/200, Y/200, pch=21, col="purple", lwd=0.6, bg=8, cex=0.8)
#dev.off()

#pdf("../figures/RVMNOTUSEDimage.pdf", useDingbats=FALSE)
plot(X/200, Y/200, type="n",
     xlab="FRACTION EASTING RANGE", ylab="FRACTION NORTHING RANGE")
image(M, add=TRUE, col=terrain.colors(20))
points(X/200, Y/200, pch=21, col="purple", lwd=0.6, bg=8, cex=0.8)
#dev.off()

Xs <- 1:200
#pdf("../draftfigures/RVMNOTUSEDhorzmarginA.pdf", useDingbats=FALSE)
set <- Y >= 50 & Y <= 100;
Xp <- X[set]; Zp <- Z[set]
plot(Xp,Zp, type="n",
     xlab="FRACTION EASTING RANGE",
     ylab="VERTICAL DISTANCE")
gam <- gam(Zp~s(Xp))
set.seed(seed)
svm <- ksvm(Zp~Xp, epsilon=0.3); six <- SVindex(svm)
points(Xp[-six], Zp[-six], col=1, lwd=0.7, pch=1,  cex=0.7)
points(Xp[ six], Zp[ six], col=grey(0.4),  pch=16, cex=0.7)
lines(Xs, predict(gam, data.frame(Xp=Xs)), col=4, lwd=2)
lines(Xs, predict(svm, data.frame(Xp=Xs)), col=2, lwd=2)
#dev.off()

#pdf("../draftfigures/RVMNOTUSEDhorzmarginB.pdf", useDingbats=FALSE)
skip <- Xp >= 100 & Xp <= 160
nX <- X[! skip]; nZ <- Zp[! skip]
plot(nX,nZ, type="n",
     xlab="FRACTION EASTING RANGE",
     ylab="VERTICAL DISTANCE")
set.seed(seed)
gam <- gam(nZ~s(nX))
set.seed(seed)
svm <- ksvm(nZ~nX, epsilon=0.3); six <- SVindex(svm)
points(nX[-six], nZ[-six], col=1, lwd=0.7, pch=1,  cex=0.7)
points(nX[ six], nZ[ six], col=grey(0.4),  pch=16, cex=0.7)
lines(Xs, predict(gam, data.frame(nX=Xs)), col=4, lwd=2)
lines(Xs, predict(svm, data.frame(nX=Xs)), col=2, lwd=2)
#dev.off()

# GAM1d <- gam(Z~s(X)+s(Y)); plot(GAM1d, residuals=TRUE)

#pdf("../draftfigures/RVMNOTUSEDgam2d.pdf", useDingbats=TRUE)
set.seed(seed)
GAM2d <- gam(Z~s(X,Y), data=data.frame(X=X, Y=Y))
NSE.gam2d.whole  <- NSE(Z,predict(GAM2d))
RMSE.gam2d.whole <- sqrt(mean(( predict(GAM2d)- Z)^2))
message("gam2d.whole: ", round(NSE.gam2d.whole, digits=3), ", ",
                         round(RMSE.gam2d.whole, digits=3))
par(las=1)
plot(GAM2d, xlab="", ylab="")
mtext("EASTING, LENGTH UNITS",  side=1, line=3); par(las=0)
mtext("NORTHING, LENGTH UNITS", side=2, line=3); par(las=1)
set.seed(seed)
SVM <- ksvm(Z~X+Y, epsilon=0.1)
pred.z <- predict(SVM, data.frame(X=X, Y=Y))
six <- SVindex(SVM); svm.n <- length(six)
points(X,       Y,       col=8, pch=16, cex=0.6)
points(X[ six], Y[ six], col=2, cex=0.8)
NSE.svm.whole  <- NSE(Z,pred.z)
RMSE.svm.whole <- sqrt(mean(( pred.z - Z)^2))
message("SVM Whole (seed, NSE, RMSE, n): ", seed, ", ",
         round(NSE.svm.whole, digits=3), ", ",
         round(RMSE.svm.whole, digits=3), ", ", svm.n, " {tab:wholemodel}")

set.seed(seed)
RVM <- rvm(Z~X+Y, kpar=list(sigma=0.01))
rix <- RVindex(RVM); rvm.n <- length(rix)
points(X[rix], Y[rix], col=4, cex=1.2)
pred.z <- predict(RVM, data.frame(X=X, Y=Y))
NSE.rvm.whole  <- NSE(Z,pred.z)
RMSE.rvm.whole <- sqrt(mean(( pred.z - Z)^2))
message("RVM Whole (seed, NSE, RMSE, n): ", seed, ", ",
         round(NSE.rvm.whole, digits=3), ", ",
         round(RMSE.rvm.whole, digits=3), ", ", rvm.n, " {tab:wholemodel}")
#dev.off()


#pdf("../draftfigures/RVMNOTUSEDsvmrvmgam.pdf", useDingbats=TRUE)
opts <- par(no.readonly = TRUE)
layout(matrix(c(1,2), nrow=2, ncol=1))
par(mar = c(4,4,.5,1))
plot(predict(GAM2d), predict(SVM, data.frame(X=X, Y=Y)),
     col=2, lwd=0.8,
     xlab="GAM 2D PREDICTION", ylab="SVM PREDICTION")
abline(0,1); text(20,120, "(A)")
plot(predict(GAM2d), predict(RVM, data.frame(X=X, Y=Y)),
     col=4, lwd=0.8,
     xlab="GAM 2D PREDICTION", ylab="RVM PREDICTION")
abline(0,1); text(20,120, "(B)")
par(opts)
#dev.off()

#save(Z,X,Y,M,nsim, file="ProbSVMexp.RData")


set.seed(seed)
GAM2d.lite <- gam(Z[six]~s(X,Y), data=data.frame(X=X[ six], Y=Y[ six]))
GAM2d.RMSE.six.mod <- sqrt(mean((predict(GAM2d.lite) - Z[six])^2))
GAM2d.NSE.six.mod  <- NSE(Z[six],predict(GAM2d.lite))
message("GAM2d.six.mod (seed, NSE, RMSE, n): ", seed, ", ",
         round(GAM2d.NSE.six.mod, digits=3), ", ",
         round(GAM2d.RMSE.six.mod, digits=3), ", ", length(Z[six]), " {tab:submodel}")
pred.z <- predict(GAM2d.lite, newdata=data.frame(X=X[-six], Y=Y[-six]))
GAM2d.RMSE.six.out <- sqrt(mean((pred.z - Z[-six])^2))
GAM2d.NSE.six.out  <- NSE(Z[-six],pred.z)
message("GAM2d.six.out (seed, NSE, RMSE, n): ", seed, ", ",
         round(GAM2d.NSE.six.out, digits=3), ", ",
         round(GAM2d.RMSE.six.out, digits=3), ", ", length(Z[-six]), " {tab:submodel}")


all <- sort(unique(c(rix,six)))
GAM2d.lite <- gam(Z[all]~s(X,Y), data=data.frame(X=X[all], Y=Y[all]))
GAM2d.RMSE.all.mod <- sqrt(mean((predict(GAM2d.lite) - Z[all])^2))
GAM2d.NSE.all.mod  <- NSE(Z[all],predict(GAM2d.lite))
message("GAM2d.all.out (seed, NSE, RMSE, n): ", seed, ", ",
        round(GAM2d.NSE.all.mod, digits=3), ", ",
        round(GAM2d.RMSE.all.mod, digits=3), ", ", length(Z[all]), " {tab:submodel}")
pred.z <- predict(GAM2d.lite, newdata=data.frame(X=X[-all], Y=Y[-all]))
GAM2d.RMSE.all.out <- sqrt(mean((pred.z- Z[-all])^2))
GAM2d.NSE.all.out  <- NSE(Z[-all],pred.z)
message("GAM2d.all.out (seed, NSE, RMSE, n): ", seed, ", ",
         round(GAM2d.NSE.all.out, digits=3), ", ",
         round(GAM2d.RMSE.all.out, digits=3), ", ", length(Z[-all]), " {tab:submodel}")


set.seed(seed)
RAND.RMSE.six <- RAND.NSE.six <- rep(NA, svm.n)
for(i in 1:nsim) {
  if(as.logical(length(grep("00$", i)))) message(i,"-", appendLF=FALSE)
  set <- sample(1:nsim, size=svm.n, replace=FALSE)
  gam.rand <- gam(Z~s(X,Y), data=data.frame(X=X[set], Y=Y[set], Z=Z[set]))
  pred.z <- predict(gam.rand, newdata=data.frame(X=X[-set], Y=Y[-set]))
  RAND.RMSE.six[i] <- sqrt(mean((pred.z - Z[-set])^2))
  RAND.NSE.six[i] <- NSE(Z[-set],pred.z)
}
message("done")

message("RAND.NSE.six:  ", seed, ", n=", svm.n, ", ",
         paste(round(summary(RAND.NSE.six),  digits=3), collapse=", "),
        " {tab:simrmsense}")
message("RAND.RMSE.six: ", seed, ", n=", svm.n, ", ",
         paste(round(summary(RAND.RMSE.six), digits=3), collapse=", "),
        " {tab:simrmsense}")


set.seed(seed)
RAND.RMSE.all <-  RAND.NSE.all <- rep(NA, length(all))
for(i in 1:nsim) {
  if(as.logical(length(grep("00$", i)))) message(i,"-", appendLF=FALSE)
  set <- sample(1:nsim, size=length(all), replace=FALSE)
  gam.rand <- gam(Z~s(X,Y), data=data.frame(X=X[set], Y=Y[set], Z=Z[set]))
  pred.z <- predict(gam.rand, newdata=data.frame(X=X[-set], Y=Y[-set]))
  RAND.RMSE.all[i] <- sqrt(mean((pred.z - Z[-set])^2))
  RAND.NSE.all[i] <- NSE(Z[-set],pred.z)
}
message("done")
message("RAND.NSE.all:  ", seed, ", n=", length(all), ", ",
        paste(round(summary(RAND.NSE.all),  digits=3), collapse=", "),
        " {tab:simrmsense}")
message("RAND.RMSE.all: ", seed, ", n=", length(all), ", ",
        paste(round(summary(RAND.RMSE.all), digits=3), collapse=", "),
        " {tab:simrmsense}")


stop()



set.seed(seed)
first <- TRUE
for(i in 1:nsim) {
  if(as.logical(length(grep("00$", i)))) message(i, "-", appendLF=FALSE)
  SVM <- ksvm(Z[-i]~X[-i]+Y[-i])
  six <- SVindex(SVM)
  if(first) { first <- FALSE
    svm <- data.frame(index=six)
  } else {
    svm <- rbind(svm, data.frame(index=six))
  }
}
message("done")
svm <- aggregate(svm, by=list(svm$index), length)
svm$count <- svm$index; svm$index <- svm$Group.1; svm$Group.1 <- NULL
plot(qnorm(lmomco::pp(svm$count)), sort(svm$count), log="y")


set.seed(seed)
first <- TRUE
for(i in 1:nsim) {
  if(as.logical(length(grep("00$", i)))) message(i, "-", appendLF=FALSE)
  RVM <- rvm(Z[-i]~X[-i]+Y[-i])
  rix <- RVindex(RVM)
  if(first) { first <- FALSE
    rvm <- data.frame(index=rix)
  } else {
    rvm <- rbind(rvm, data.frame(index=rix))
  }
}
message("done")
rvm <- aggregate(rvm, by=list(rvm$index), length)
rvm$count <- rvm$index; rvm$index <- rvm$Group.1; rvm$Group.1 <- NULL
plot(qnorm(lmomco::pp(rvm$count)), sort(rvm$count), log="y")


plot(GAMd2, scheme=3)
points(X, Y, col=8, pch=16, cex=0.6)
points(X[svm$index[svm$count <=    1]],
       Y[svm$index[svm$count <=    1]], col=6, pch=16, cex=1.5)
points(X[svm$index[svm$count == 800]],
       Y[svm$index[svm$count == 800]], col=5, pch=16, cex=1.5)

set <- svm$index[svm$count > 800]
gam.rand <- gam(Z~s(X,Y), data=data.frame(X=X[set], Y=Y[set], Z=Z[set]))
pred.z <- predict(gam.lite, newdata=data.frame(X=X[-set], Y=Y[-set]))
RAND.RMSE.every <- sqrt(mean((pred.z - Z[-set])^2))
RAND.NSE.every  <- NSE(Z[-set],pred.z)



grid <- 1:200
P.mod <- matrix(nrow=length(grid), ncol=length(grid))
for(i in 1:length(grid)) {
 for(j in 1:length(grid)) {
  P.mod[i,j] <- predict(GAM2d, data.frame(X=grid[i], Y=grid[j]))
 }
}
grid <- 1:200
P.all <- matrix(nrow=length(grid), ncol=length(grid))
for(i in 1:length(grid)) {
 for(j in 1:length(grid)) {
  P.all[i,j] <- predict(GAM2d.lite, data.frame(X=grid[i], Y=grid[j]))
 }
}

#pdf("../draftfigures/RVMNOTUSEDsurface_modminusall.pdf", useDingbats=FALSE)
plot(X/200, Y/200, type="n",
     xlab="FRACTION EASTING RANGE", ylab="FRACTION NORTHING RANGE")
contour(P.mod-P.all, lwd=1.5,labcex = 0.7, nlevels=8, add=TRUE, axes=FALSE)
points(X/200, Y/200, pch=21, col="purple", lwd=0.6, bg=8, cex=0.8)
#dev.off()

#pdf("../draftfigures/RVMNOTUSEDcontour.pdf", useDingbats=FALSE)
contour(P.mod-P.all)
#dev.off()





