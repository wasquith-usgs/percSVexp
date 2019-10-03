### ----------------------------------------------------------------
library(kernlab)
library(mgcv)
seed <- 1 # notice that set.seed(62) is used about 25 lines later, here this is a reminder to the reader

# STEP 1.
nsim <- 20 # set the number of simulations

NSE <- function(obs,sim) { # Nash-Sutcliffe Efficiency
   1 - (sum((obs-sim)^2)/sum((obs-mean(obs))^2))
}

xymod <- function(x,y, datum=70) {
  i <- 20*sin(+1.75*pi*(x-y)/100)
  j <- 20*cos(-1.50*pi*   y /100)
  z <- i + j
  err <- rnorm(length(x), mean=0, sd=2)
  drift <- 0.001*x^2 - 0.1*y
  z <- z + drift + err + datum
  return(z)
}

# STEP 2.
grid <- 1:200
M <- matrix(nrow=length(grid), ncol=length(grid))
for(i in 1:length(grid)) {
  for(j in 1:length(grid)) {
    M[i,j] <-  xymod(grid[i], grid[j])
  }
}

# STEP 3.
set.seed(62)
nobs <- 800
X <- runif(nobs, min=0, 200)
Y <- runif(nobs, min=0, 200)
Z <- xymod(X,Y)
set.seed(seed)

# STEP 4.
pdf("../draftfigures/fig01_rawsurface.pdf", useDingbats=FALSE)
  opts <- par(no.readonly = TRUE); par(las=1)
  plot(X/200, Y/200, type="n",
       xlab="FRACTION EASTING RANGE", ylab="FRACTION NORTHING RANGE")
  polygon(c(0,1,1,0,0), c(50,50,100,100,50)/200, border=2, col=8, lwd=2, lty=2)
  points(X/200, Y/200, pch=21, col="purple", lwd=0.6, bg=8, cex=0.8)
  contour(M, lwd=1.5, labcex = 1.2, nlevels=8, add=TRUE, axes=FALSE, col=grey(0.3))
  par(opts)
dev.off()

pdf("../draftfigures/fig02_rawlegend.pdf", useDingbats=FALSE)
  opts <- par(no.readonly = TRUE); par(las=1)
  cols <- terrain.colors(20)
  quas <- quantile(Z, probs=(1:20)/20)
  plot(rep(1,20), 1:20, col=rev(cols), pch=15, cex=3.2,
       xaxt="n", yaxt="n", xlab="", ylab="")
  for(i in 1:20) {
    text(1.1, i, round(quas[i], digits=0))
  }
  mtext("Approximate break points\nin the Z direction")
  par(opts)
dev.off()

pdf("../draftfigures/fig02_rawimage.pdf", useDingbats=FALSE)
  opts <- par(no.readonly = TRUE); par(las=1)
  plot(X/200, Y/200, type="n",
       xlab="FRACTION EASTING RANGE", ylab="FRACTION NORTHING RANGE")
  image(M, add=TRUE, col=terrain.colors(20))
  polygon(c(0,1,1,0,0), c(50,50,100,100,50)/200, border=2, col=NA, lwd=2, lty=2)
  points(X/200, Y/200, pch=21, col="purple", lwd=0.6, bg=8, cex=0.8)
  par(opts)
dev.off()
### ----------------------------------------------------------------




# STEP 5
Xs <- 1:200
pdf("../draftfigures/fig03_horzmarginA.pdf", useDingbats=FALSE)
  opts <- par(no.readonly = TRUE); par(las=1)
  notskip <- Y >= 50 & Y <= 100;
  Xp <- X[notskip]; Zp <- Z[notskip]
  plot(Xp,Zp, type="n",
       xlab="EASTING DISTANCE",
       ylab="VERTICAL DISTANCE")
  lines(par()$usr[1:2], rep(mean(Zp), 2), lty=2)
  gam <- mgcv::gam(Zp~s(Xp, bs="tp"))
  set.seed(seed)
  svm <- kernlab::ksvm(Zp~Xp, epsilon=0.3); six <- kernlab::SVindex(svm)
  lines(Xs, predict(svm, data.frame(Xp=Xs)), col=2, lwd=2)
  points(Xp[-six], Zp[-six], col=1, lwd=0.7, pch=1,  cex=1.1)
  points(Xp[ six], Zp[ six], col=grey(0.4),  pch=16, cex=0.8)
  lines(Xs, predict(gam, data.frame(Xp=Xs)), col=4, lwd=2)
  par(opts)
dev.off()
message("Global vertical mean in partition: ", round(mean(Zp), digits=3))


# STEP 6
pdf("../draftfigures/fig04_horzmarginB.pdf", useDingbats=FALSE)
  opts <- par(no.readonly = TRUE); par(las=1)
  skip <- Xp >= 100 & Xp <= 160
  nX <- Xp[! skip]; nZ <- Zp[! skip]
  plot(nX,nZ, type="n",
       xlab="EASTING DISTANCE",
       ylab="VERTICAL DISTANCE")
  lines(par()$usr[1:2], rep(mean(nZ), 2), lty=2)
  set.seed(seed)
  gam <- mgcv::gam(nZ~s(nX, bs="tp"))
  set.seed(seed)
  svm <- kernlab::ksvm(nZ~nX, epsilon=0.3); six <- kernlab::SVindex(svm)
  lines(Xs, predict(svm, data.frame(nX=Xs)), col=2, lwd=2)
  points(nX[-six], nZ[-six], col=1, lwd=0.7, pch=1,  cex=1.1)
  points(nX[ six], nZ[ six], col=grey(0.4),  pch=16, cex=0.8)
  lines(Xs, predict(gam, data.frame(nX=Xs)), col=4, lwd=2)
  par(opts)
dev.off()
message("Global vertical mean in cut partition: ", round(mean(nZ), digits=3))

# STEP 7
txt <- "global mean of vertical distance"
pdf("../draftfigures/fig05_horzmarginC.pdf", useDingbats=FALSE)
  opts <- par(no.readonly = TRUE); par(las=1)
  skip <- Xp >= 100 & Xp <= 160
  nX <- Xp[! skip]; nZ <- Zp[! skip]
  plot(nX,nZ, type="n",
       xlab="EASTING DISTANCE",
       ylab="VERTICAL DISTANCE")
  lines(par()$usr[1:2], rep(mean(nZ), 2), lty=2)
  points(nX, nZ, col=1,  pch=16, cex=0.9)
  set.seed(62)
  for(i in 1:300) {
    svm <- kernlab::ksvm(nZ~nX, epsilon=0.3); six <- kernlab::SVindex(svm)
    lines(Xs, predict(svm, data.frame(nX=Xs)), col=rgb(1,0,0,.1), lwd=1.1)
  }
  set.seed(seed)
  svm <- kernlab::ksvm(nZ~nX, epsilon=0.3); six <- kernlab::SVindex(svm)
  lines(Xs, predict(svm, data.frame(nX=Xs)), col=rgb(0.1,0.8,.1), lwd=3)
  text(135, 40, txt, cex=0.9)
  arrows(100, 42, 100, mean(nZ), angle=15)
  par(opts)
dev.off()


# STEP 8
# GAM1d <- gam(Z~s(X)+s(Y)); plot(GAM1d, residuals=TRUE) # just an experiment on this line
pdf("../draftfigures/fig06_gam2d.pdf", useDingbats=FALSE)
  opts <- par(no.readonly = TRUE); par(las=1)
  set.seed(seed)
  GAM2d <- mgcv::gam(Z~s(X,Y, bs="tp"), data=data.frame(X=X, Y=Y))
  message("gam2d.whole: gam intercept is ",round(GAM2d$coefficients[1], digits=2))
  NSE.gam2d.whole  <- NSE(Z,predict(GAM2d))
  RMSE.gam2d.whole <- sqrt(mean(( predict(GAM2d)- Z)^2))
  message("gam2d.whole: ", round(NSE.gam2d.whole, digits=3), ", ",
                           round(RMSE.gam2d.whole, digits=3), " {tab:wholemodel}")
  mgcv::plot.gam(GAM2d, xlab="", ylab="")
  mtext("EASTING DISTANCE",  side=1, line=3); par(las=0)
  mtext("NORTHING DISTANCE", side=2, line=3); par(las=1)
  set.seed(seed)
  SVM <- kernlab::ksvm(Z~X+Y, epsilon=0.1)
  pred.z <- predict(SVM, data.frame(X=X, Y=Y))
  six <- kernlab::SVindex(SVM); svm.n <- length(six)
  points(X,       Y,       col=8, pch=16, cex=0.6)
  points(X[ six], Y[ six], col=2, cex=0.8)
  NSE.svm.whole  <- NSE(Z,pred.z)
  RMSE.svm.whole <- sqrt(mean(( pred.z - Z)^2))
  message("SVM Whole (seed, NSE, RMSE, n): ", seed, ", ",
           round(NSE.svm.whole,  digits=3), ", ",
           round(RMSE.svm.whole, digits=3), ", ", svm.n, " {tab:wholemodel}")
  par(opts)
dev.off()


pdf("../draftfigures/fig07_svmgam.pdf", useDingbats=FALSE)
  opts <- par(no.readonly = TRUE); par(las=1)
  plot(predict(GAM2d), predict(SVM, data.frame(X=X, Y=Y)),
       col="#66A48B", lwd=0.8,
       xlab="GAM PREDICTION, UNITLESS", ylab="SVM PREDICTION, UNITLESS")
  abline(0,1)
  par(opts)
dev.off()

save(Z,X,Y,M, file="ProbSVMexp.RData")




# STEP 9
set.seed(seed)
GAM2d.lite <- mgcv::gam(Z[six]~s(X,Y, bs="tp"), data=data.frame(X=X[ six], Y=Y[ six]))
GAM2d.RMSE.six.mod <- sqrt(mean((predict(GAM2d.lite) - Z[six])^2))
GAM2d.NSE.six.mod  <- NSE(Z[six],predict(GAM2d.lite))
#message("GAM2d.six.mod (seed, NSE, RMSE, n): ", seed, ", ",
#         round(GAM2d.NSE.six.mod, digits=3), ", ",
#         round(GAM2d.RMSE.six.mod, digits=3), ", ", length(Z[six]), " {tab:submodel}")
pred.z <- predict(GAM2d.lite, newdata=data.frame(X=X[-six], Y=Y[-six]))
GAM2d.RMSE.six.out <- sqrt(mean((pred.z - Z[-six])^2))
GAM2d.NSE.six.out  <- NSE(Z[-six],pred.z)
#message("GAM2d.six.out (seed, NSE, RMSE, n): ", seed, ", ",
#         round(GAM2d.NSE.six.out, digits=3), ", ",
#         round(GAM2d.RMSE.six.out, digits=3), ", ", length(Z[-six]), " {tab:submodel}")

message("Now computing grid level whole model results")

P.lite <- matrix(nrow=length(grid), ncol=length(grid))
for(k in 1:length(grid)) {
  P.lite[k,]  <- predict(GAM2d.lite,
                         data.frame(X=rep(grid[k],200), Y=grid))
}
GAM2d.invlite <- mgcv::gam(Z[-six]~s(X,Y, bs="tp"), data=data.frame(X=X[-six], Y=Y[-six]))

P.invlite <- matrix(nrow=length(grid), ncol=length(grid))
for(k in 1:length(grid)) {
  P.invlite[k,] <- predict(GAM2d.invlite,
                            data.frame(X=rep(grid[k],200), Y=grid))
}
message("Grid Level +n^svm Model (NSE, RMSE): n=",svm.n, ", ",
        round(NSE(M,P.lite), digits=3), ", ",
        round(sqrt(mean((M-P.lite)^2)), digits=3),
        " {tab:submodel}")
message("Grid Level -n^svm Model NSE: n=",nobs-svm.n, ", ",
        round(NSE(M,P.invlite), digits=3), ", ",
        round(sqrt(mean((M-P.invlite)^2)), digits=3),
        " {tab:submodel}")



# STEP 10
set.seed(seed)
RAND.RMSE.rand <- RAND.NSE.rand <- rep(NA, nsim); ix <- 1:800
for(j in 1:nsim) {
  if(as.logical(length(grep("0$", j)))) message(j,"-", appendLF=FALSE)
  set <- sample(ix, size=svm.n, replace=FALSE)
  gam.rand <- mgcv::gam(Z~s(X,Y, bs="tp"), data=data.frame(X=X[set], Y=Y[set], Z=Z[set]))
  P.lite <- matrix(nrow=length(grid), ncol=length(grid))
  for(k in 1:length(grid)) {
    P.lite[k,] <- predict(gam.rand,
                          data.frame(X=rep(grid[k],200), Y=grid))
  }
  #pred.z <- predict(gam.rand, newdata=data.frame(X=X[-set], Y=Y[-set]))
  RAND.RMSE.rand[j] <- sqrt(mean((M-P.lite)^2))#sqrt(mean((pred.z - Z[-set])^2))
  RAND.NSE.rand[j]  <- NSE(M,P.lite)#NSE(Z[-set],pred.z)
}
message("done")
message("RAND.NSE.rand:  ", seed, ", n=", svm.n, ", ",
         paste(round(summary(RAND.NSE.rand),  digits=3), collapse=", "),
        " {tab:simrmsense}")
message("RAND.RMSE.rand: ", seed, ", n=", svm.n, ", ",
         paste(round(summary(RAND.RMSE.rand), digits=3), collapse=", "),
        " {tab:simrmsense}")


set.seed(seed)
RAND.RMSE.six <- RAND.NSE.six <- SVM.Ns <- rep(NA, nsim)
for(j in 1:nsim) {
  if(as.logical(length(grep("0$", j)))) message(j,"-", appendLF=FALSE)
  svm <- kernlab::ksvm(Z~X+Y, epsilon=0.1); set <- kernlab::SVindex(svm)
  SVM.Ns[j] <- length(set)
  gam.svm <- mgcv::gam(Z~s(X,Y, bs="tp"), data=data.frame(X=X[set], Y=Y[set], Z=Z[set]))
  P.lite <- matrix(nrow=length(grid), ncol=length(grid))
  for(k in 1:length(grid)) {
    P.lite[k,] <- predict(gam.svm,
                          data.frame(X=rep(grid[k],200), Y=grid))
  }
  #pred.z <- predict(gam.rand, newdata=data.frame(X=X[-set], Y=Y[-set]))
  RAND.RMSE.six[j] <- sqrt(mean((M-P.lite)^2))#sqrt(mean((pred.z - Z[-set])^2))
  RAND.NSE.six[j]  <- NSE(M,P.lite)#NSE(Z[-set],pred.z)
}
message("done")
message("RAND.NSE.six:  ", seed, ", n=", svm.n, ", ",
         paste(round(summary(RAND.NSE.six),  digits=3), collapse=", "),
        " {tab:simrmsense}")
message("RAND.RMSE.six: ", seed, ", n=", svm.n, ", ",
         paste(round(summary(RAND.RMSE.six), digits=3), collapse=", "),
        " {tab:simrmsense}")
message("SVM.Ns with simulations:",
         paste(round(summary(SVM.Ns), digits=3), collapse=", "),
        "{tab:distsamplesize}")

#stop()



set.seed(seed)
first <- TRUE; ix <- 1:length(Z) # seq of overall indices
message("Starting simulation chucks totaling ",nsim," and this takes time!")
for(j in 1:nsim) {
  message("simulation chuck: ",j)
  for(i in ix) {
    if(as.logical(length(grep("00$", i)))) message(i, "-", appendLF=FALSE)
    SVM <- kernlab::ksvm(Z[-i]~X[-i]+Y[-i]) # note the subsettings of the ix to
    tix <- ix[-i]; six <- tix[kernlab::SVindex(SVM)] # leave-one-out and then six'ing
    if(first) { first <- FALSE
      svm <- data.frame(index=six)
    } else {
      svm <- rbind(svm, data.frame(index=six))
    }
  }
  message("done")
}
message("done")
svm <- aggregate(svm, by=list(svm$index), length)
svm$svm_count <- svm$index; svm$index <- svm$Group.1; svm$Group.1 <- NULL
svm$svm_ratio <- svm$svm_count/max(svm$svm_count)

svmck <- rbind(data.frame(index=svm$index), data.frame(index=ix))
svmck <- aggregate(svmck, by=list(svmck$index), length)
indices_never_used <- svmck$Group.1[svmck$index == 1]
svm2 <- rbind(svm, data.frame(index=indices_never_used, svm_count=0, svm_ratio=0))
svm2 <- svm2[order(svm2$index),]
svm <- svm2

txt <- "Color hue is prorated from\nred to blue based on\nnonexceedance probability."
pdf("../draftfigures/fig08_svmvecpp.pdf", useDingbats=FALSE)
  opts <- par(no.readonly = TRUE); par(las=1)
  tmp <- svm[order(svm$svm_ratio),]
  plot(lmomco::pp(tmp$svm_ratio, sort=FALSE),
                  tmp$svm_ratio*100,
       xlab="NONEXCEEDANCE PROBABILITY", lwd=0.8, type="p",
       ylab="SUPPORT VECTOR PERCENTAGE", xlim=c(0,1),
       col=rgb(1-tmp$svm_ratio,0,tmp$svm_ratio))
  text(0.3, 70, txt, cex=0.9)
  par(opts)
dev.off()


plot(GAM2d, scheme=3)
points(X, Y, col=8, pch=16, cex=0.6)
points(X[svm$index[svm$svm_count <=   5]],
       Y[svm$index[svm$svm_count <=   5]], col=6, pch=16, cex=1.5)
points(X[svm$index[svm$svm_count >= 795]],
       Y[svm$index[svm$svm_count >= 795]], col=5, pch=16, cex=1.5)

set <- svm$index[svm$svm_count > 500]
gam.lite <- mgcv::gam(Z~s(X,Y, bs="tp"), data=data.frame(X=X[set], Y=Y[set], Z=Z[set]))
pred.z <- predict(gam.lite, newdata=data.frame(X=X[-set], Y=Y[-set]))
RAND.RMSE.every <- sqrt(mean((pred.z - Z[-set])^2))
RAND.NSE.every  <- NSE(Z[-set],pred.z)


P.mod <- matrix(nrow=length(grid), ncol=length(grid))
for(k in 1:length(grid)) {
  P.mod[k,]  <- predict(GAM2d,
                       data.frame(X=rep(grid[k],200), Y=grid))
}
P.six <- matrix(nrow=length(grid), ncol=length(grid))
for(k in 1:length(grid)) {
 P.six[k,]  <- predict(GAM2d.lite,
                        data.frame(X=rep(grid[k],200), Y=grid))
}

# figure not used in the paper
#pdf("../draftfigures/NOTUSEDsurface_modminusall.pdf", useDingbats=FALSE)
  opts <- par(no.readonly = TRUE); par(las=1)
  plot(X/200, Y/200, type="n",
     xlab="FRACTION EASTING RANGE", ylab="FRACTION NORTHING RANGE")
  contour(P.mod-P.six, lwd=1.5,labcex = 0.7, nlevels=8, add=TRUE, axes=FALSE)
  points(X/200, Y/200, pch=21, col="purple", lwd=0.6, bg=8, cex=0.8)
  par(opts)
#dev.off()

#pdf("../draftfigures/NOTUSEDcontour.pdf", useDingbats=FALSE)
contour(P.mod-P.six) # figure not used in the paper
#dev.off()




