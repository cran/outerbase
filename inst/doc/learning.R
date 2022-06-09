## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(outerbase)

## -----------------------------------------------------------------------------
om = new(outermod)
d = 3
setcovfs(om, rep("mat25pow",d))
knotlist = list(seq(0.01,0.99,by=0.02),
                seq(0.01,0.99,by=0.02),
                seq(0.01,0.99,by=0.02))
setknot(om, knotlist)

## -----------------------------------------------------------------------------
sampsize = 30
design1d = seq(1/(2*sampsize),1-1/(2*sampsize),1/sampsize)
x = cbind(design1d,sample(design1d),sample(design1d))
ob = new(outerbase, om, x)
basis_func0 = ob$getbase(1)
matplot(x[,1],basis_func0[,1:4], 
        type='l', ylab="func", xlab="first dim")

## -----------------------------------------------------------------------------
hyp0 = gethyp(om)
hyp0[2] = 3 #changing the power on first parameter
om$updatehyp(hyp0)
ob$build() #rebuild after updatehyp

## ----fig.show="hold", out.width="45%", fig.width=4, fig.height=4--------------
basis_func1 = ob$getbase(1)
matplot(x[,1],basis_func0[,1:4], 
        type='l', ylab="func", xlab="first dim",
        main="original hyperparameters")
matplot(x[,1],basis_func1[,1:4], 
        type='l', ylab="func", xlab="first dim",
        main="new hyperparameters")

## -----------------------------------------------------------------------------
y = obtest_borehole3d(x)

## -----------------------------------------------------------------------------
gethyp(om)
hyp0 = c(-0.5,0,-0.5,0,-0.5,0)
om$updatehyp(hyp0)

## -----------------------------------------------------------------------------
terms = om$selectterms(60)

## -----------------------------------------------------------------------------
loglik = new(loglik_std, om, terms, y, x) 
coeff0 = rep(0,loglik$nterms)
loglik$update(coeff0) # update it to get gradients
loglik$val
head(loglik$grad) # dim 60 for number of coeffients

## -----------------------------------------------------------------------------
logpr = new(logpr_gauss, om, terms)

## -----------------------------------------------------------------------------
logpdf = new(lpdfvec, loglik, logpr)
para0 = getpara(logpdf)
para0
para0[2] = 4
logpdf$updatepara(para0)
getpara(logpdf)

## -----------------------------------------------------------------------------
logpdf$optnewton()

## -----------------------------------------------------------------------------
testsampsize = 1000
xtest = matrix(runif(testsampsize*d),ncol=d)
ytest = obtest_borehole3d(xtest)

## ----fig.show="hold", out.width="45%", fig.width=4, fig.height=4--------------
predt = new(predictor,loglik)
predt$update(xtest)
yhat = as.vector(predt$mean())
varpred = as.vector(predt$var())

plot(yhat,ytest, xlab="prediction", ylab="actual")
hist((ytest-yhat)/sqrt(varpred),
     main="standarized test residuals",
     xlab = "standarized test residuals")

## -----------------------------------------------------------------------------
logpdf$optnewton()
logpdf$gradhyp    # dim 6 for all hyperparameter
logpdf$gradpara   # dim 2 since 2 parameters

## -----------------------------------------------------------------------------
totobj = function(parlist) { #my optimization function for tuning
  regpara = logpdf$paralpdf(parlist$para) # get regularization for lpdf
  reghyp = om$hyplpdf(parlist$hyp) # get regularization for om
  if(is.finite(regpara) && is.finite(reghyp)) { # if they pass
    om$updatehyp(parlist$hyp)        # update hyperparameters
    logpdf$updateom()             # update the outerbase inside
    logpdf$updatepara(parlist$para)  # update parameter
    logpdf$optnewton()            # do opt
    
    gval = parlist #match structure
    gval$hyp = -logpdf$gradhyp-om$hyplpdf_grad(parlist$hyp)
    gval$para = -logpdf$gradpara-logpdf$paralpdf_grad(parlist$para)
    list(val = -logpdf$val-reghyp-regpara, gval = gval)
  } else list(val = Inf, gval = parlist) }

## -----------------------------------------------------------------------------
parlist = list(para = getpara(logpdf), hyp = gethyp(om))
totobj(parlist)
opth = BFGS_std(totobj, parlist, verbose=3) #

## -----------------------------------------------------------------------------
totobj(opth$parlist)

## -----------------------------------------------------------------------------
opth = BFGS_lpdf(om, logpdf, 
                 parlist=parlist, 
                 verbose = 3, newt= TRUE)  

## ----fig.show="hold", out.width="45%", fig.width=4, fig.height=4--------------
predtt = new(predictor,loglik)
predtt$update(xtest)
yhat = as.vector(predtt$mean())
varpred = as.vector(predtt$var())

plot(ytest,yhat)
hist((ytest-yhat)/sqrt(varpred), main="standarized test residuals",
     xlab = "standarized test residuals")

