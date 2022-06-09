## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(outerbase)

## -----------------------------------------------------------------------------
sampsize = 500
d = 8
x = matrix(runif(sampsize*d),ncol=d)
y = obtest_borehole8d(x)

## -----------------------------------------------------------------------------
om = new(outermod)
setcovfs(om, rep("mat25pow",8))
knotlist = list();
for(k in 1:d) knotlist[[k]] = seq(0.01,1,by=0.025)
setknot(om, knotlist) #40 knot point for each dim

## -----------------------------------------------------------------------------
p = 250
terms = om$selectterms(p)

## -----------------------------------------------------------------------------
loglik_slow = new(loglik_std, om, terms, y, x) 
logpr_slow = new(logpr_gauss, om, terms)
logpdf_slow = new(lpdfvec, loglik_slow, logpr_slow)

## -----------------------------------------------------------------------------
logpdf_slow$optnewton()

## -----------------------------------------------------------------------------
loglik_fast = new(loglik_gauss, om, terms, y, x) 
logpr_fast = new(logpr_gauss, om, terms)
logpdf_fast = new(lpdfvec, loglik_fast, logpr_fast)

## ---- error=TRUE--------------------------------------------------------------
logpdf_fast$optnewton()

## -----------------------------------------------------------------------------
logpdf_fast$optcg(0.001,  # tolerance
                  100)    # max epochs

## -----------------------------------------------------------------------------
ob = new(outerbase, om, x) 
ob$nthreads

## -----------------------------------------------------------------------------
logpdf_slow$setnthreads(4)
logpdf_fast$setnthreads(4)

## -----------------------------------------------------------------------------
parlist_slow = list(para = getpara(logpdf_slow), hyp = gethyp(om))
parlist_fast = list(para = getpara(logpdf_fast), hyp = gethyp(om))

## -----------------------------------------------------------------------------
xtest = matrix(runif(1000*d),ncol=d) #prediction points
ytest =  obtest_borehole8d(xtest)

## -----------------------------------------------------------------------------
ptm = proc.time()
opth = BFGS_lpdf(om, logpdf_slow, 
                 parlist=parlist_slow, newt=TRUE)    
t_slow = proc.time() - ptm
pred_slow = new(predictor,loglik_slow)
pred_slow$update(xtest)
yhat_slow = as.vector(pred_slow$mean())
print(t_slow)

## -----------------------------------------------------------------------------
ptm = proc.time()
opth = BFGS_lpdf(om, logpdf_fast, 
                 parlist=parlist_fast, newt=FALSE)  
t_fast = proc.time() - ptm
pred_fast = new(predictor,loglik_fast)
pred_fast$update(xtest)
yhat_fast = as.vector(pred_fast$mean())
print(t_fast)

## ----fig.show="hold", out.width="45%", fig.width=4, fig.height=4--------------
rmse_slow = sqrt(mean((ytest-yhat_slow)^2))
hist((ytest-yhat_slow), main=paste("slow method \n rmse:", 
                                    round(rmse_slow,3),
                                   ", time:",
                                   round(t_slow[3],2),'s'),
     xlab = "prediction residuals")
rmse_fast = sqrt(mean((ytest-yhat_fast)^2))
hist((ytest-yhat_fast), main=paste("fast method \n rmse =",
                                      round(rmse_fast,3),
                                   ", time:",
                                   round(t_fast[3],2),'s'), 
     xlab = "prediction residuals")

