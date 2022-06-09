## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(outerbase)

## -----------------------------------------------------------------------------
sampsize = 30
d = 3
design1d = seq(1/(2*sampsize),1-1/(2*sampsize),1/sampsize)
x = cbind(design1d,sample(design1d),sample(design1d))
y = obtest_borehole3d(x)

## -----------------------------------------------------------------------------
corf = new(covf_mat25)

## -----------------------------------------------------------------------------
xred = x[1:5,1]
print(corf$cov(xred,xred),3)

## -----------------------------------------------------------------------------
corf$hyp

## -----------------------------------------------------------------------------
corf$hyp = c(-0.5)
plot(x[,1],corf$cov(x[,1],0.5), type='l',
     ylab='correlation with 0.5', xlab='input')
corf$hyp = c(-0.25)
lines(x[,1],corf$cov(x[,1],0.5), type='l')
corf$hyp = c(0)
lines(x[,1],corf$cov(x[,1],0.5), type='l')

## -----------------------------------------------------------------------------
corf1 = new(covf_mat25)
corf2 = new(covf_mat25)
corf3 = new(covf_mat25)
corf1$hyp = c(-0.5) # just setting them all to the same 
corf2$hyp = c(-0.5) # hyperparameter for now
corf3$hyp = c(-0.5)

## -----------------------------------------------------------------------------
covftot = function(x1,x2){
  corf1$cov(x1[,1],x2[,1])*
  corf2$cov(x1[,2],x2[,2])*
  corf3$cov(x1[,3],x2[,3])
}
cormattot = covftot(x,x) #total correlation matrix

## -----------------------------------------------------------------------------
testsampsize = 1000
xtest = matrix(runif(testsampsize*d),ncol=d)

## -----------------------------------------------------------------------------
yhat = covftot(xtest,x) %*% solve(cormattot,y)

## ----fig.show="hold", out.width="45%", fig.width=4, fig.height=4--------------
ytest = obtest_borehole3d(xtest)
plot(yhat, ytest, ylab="actual", xlab="prediction") 
hist(ytest-yhat, main="test residuals",
     xlab = "test residuals")

## -----------------------------------------------------------------------------
sigma2hat = as.double(t(y)%*% solve(cormattot,y)/length(y))

varpred = sigma2hat*(covftot(xtest,xtest)-t(covftot(x,xtest))%*%
  solve(cormattot,covftot(x,xtest)))
hist((ytest-yhat)/sqrt(diag(varpred)),
     main="standarized test residuals",
     xlab = "standarized test residuals")

## -----------------------------------------------------------------------------
om = new(outermod)

## -----------------------------------------------------------------------------
setcovfs(om, rep("mat25",3))

## -----------------------------------------------------------------------------
knotlist = list(seq(0,1,by=0.025),
                seq(0,1,by=0.025),
                seq(0,1,by=0.025))
setknot(om, knotlist)

## -----------------------------------------------------------------------------
gethyp(om)
om$updatehyp(c(-0.5,-0.5,-0.5))
gethyp(om)

## -----------------------------------------------------------------------------
ob = new(outerbase, 
         om, # an outermod (reference only)
          x) # an input matrix

## -----------------------------------------------------------------------------
basis_func = ob$getbase(1)
matplot(x[,1],basis_func[,1:4], 
        type='l', ylab="func", xlab="first dim")

## -----------------------------------------------------------------------------
p = 60
terms = om$selectterms(p) # 60 by 3 matrix
head(terms)

## -----------------------------------------------------------------------------
covcoeff = as.vector(om$getvar(terms))

## -----------------------------------------------------------------------------
basismat = ob$getmat(terms)

termno = 5
basevec = ob$getbase(1)[,terms[termno,1]+1]*
  ob$getbase(2)[,terms[termno,2]+1]*
  ob$getbase(3)[,terms[termno,3]+1]

cbind(basevec[1:5],basismat[1:5,5]) # expect equal

## -----------------------------------------------------------------------------
cormatob = basismat%*%diag(covcoeff)%*%t(basismat)

print(round(cormattot[1:5,1:5],3)) # typical gp
print(round(cormatob[1:5,1:5],3)) # outerbase

## -----------------------------------------------------------------------------
noisevar = 10^(-4)
#posterior precision matrix of coefficients
postcov = solve(1/noisevar*t(basismat)%*%basismat+ 
                  1/sigma2hat*diag(1/covcoeff))
#posterior mean of coefficients
coeffest = postcov%*%(1/noisevar*t(basismat)%*%y)

## -----------------------------------------------------------------------------
obtest = new(outerbase, 
         om,     # same outermod 
          xtest) # new input matrix

basistest = obtest$getmat(terms)

## ----fig.show="hold", out.width="45%", fig.width=4, fig.height=4--------------
yhatob = basistest%*%coeffest
plot(yhat, ytest, main="typical gp",
     xlab="prediction", ylab="actual")
plot(yhatob, ytest, main = "outerbase equiv.",
     xlab="prediction", ylab="actual")

## ----fig.show="hold", out.width="45%", fig.width=4, fig.height=4--------------
hist(ytest-yhat, main="typical gp",
     xlab="test residuals")
hist(ytest-yhatob, main="outerbase equiv.",
     xlab="test residuals")

## ----fig.show="hold", out.width="45%", fig.width=4, fig.height=4--------------
varpredob = basistest%*%postcov%*%t(basistest) 
hist((ytest-yhat)/sqrt(diag(varpred)), main="typical gp",
     xlab="standarized test residuals")
hist((ytest-yhatob)/sqrt(diag(varpredob)), main="outerbase equiv.",
     xlab="standarized test residuals")

