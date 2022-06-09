## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(outerbase)

## -----------------------------------------------------------------------------
sampsize = 400
d = 8
x = matrix(runif(sampsize*d),ncol=d) #uniform samples
y = obtest_borehole8d(x) + 0.5*rnorm(sampsize)

## -----------------------------------------------------------------------------
listcov()

## ---- error=TRUE--------------------------------------------------------------
obmodel = obfit(x, y, covnames=rep("elephant",8))
obmodel = obfit(x, y[1:200], covnames=rep("mat25pow",5))
obmodel = obfit(x[1:2,], y[1:2], covnames=rep("mat25pow",8))
obmodel = obfit(x, y, numb = 2, covnames=rep("mat25pow",8))
obmodel = obfit(100*x, y, covnames=rep("mat25pow",8))
obmodel = obfit(0.001*x, y, covnames=rep("mat25pow",8))

## -----------------------------------------------------------------------------
ptm = proc.time()
obmodel = obfit(x, y, numb=300, covnames=rep("mat25pow",8),
                verbose = 3) 
print((proc.time() - ptm)[3])

## -----------------------------------------------------------------------------
ptm = proc.time()
obmodel = obfit(x, y, numb=300, covnames=rep("mat25pow",8),
                nthreads=1) #optional input
print((proc.time() - ptm)[3])

## -----------------------------------------------------------------------------
predtr = obpred(obmodel, x)
rmsetr = sqrt(mean((y-predtr$mean)^2))
plot(predtr$mean, y,
     main=paste("training \n RMSE = ", round(rmsetr,3)),
     xlab="prediction", ylab = "actual")

## -----------------------------------------------------------------------------
ytrue = obtest_borehole8d(x)
rmsetr = sqrt(mean((ytrue-predtr$mean)^2))
plot(predtr$mean, ytrue,
     main=paste("oracle \n RMSE = ", round(rmsetr,3)),
     xlab="prediction", ylab="actual")

## -----------------------------------------------------------------------------
xtest = matrix(runif(1000*d),ncol=d) #prediction points
ytest = obtest_borehole8d(xtest) + 0.5*rnorm(1000)

## -----------------------------------------------------------------------------
predtest = obpred(obmodel, xtest)

rmsetst = sqrt(mean((ytest-predtest$mean)^2))
plot(predtest$mean, ytest, 
     main=paste("testing \n RMSE = ", round(rmsetst,3)),
     xlab="prediction", ylab="actual")

## -----------------------------------------------------------------------------
hist((ytest-predtest$mean),
     main="testing \n  residuals", xlab="residuals")
hist((ytest-predtest$mean)/sqrt(predtest$var),
     main="testing \n standarized residuals",
     xlab="standarized residuals")

