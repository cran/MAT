score <-
function(ipar,resp,th,p,sigma,maxIter=30,conv=0.001,D=1.7,Fisher=T) {
#calculates MAP theta and SE estimates for a single person
ni<-nrow(ipar)
pos.ths<-th
iter<-0
converged=F

while (iter<maxIter && !converged) {
iter<-iter+1
pre.ths<-pos.ths

deriv1<-dLL(ipar,resp,pre.ths,p,sigma,D=D,Bayesian=T)

if (Fisher) deriv2<- -makeFI(ipar,pre.ths,p,sigma,D=D,add.sigma=T)
else deriv2<-makeHessian(ipar,resp,pre.ths,p,sigma,D=D,add.sigma=T)

delta<-solve(deriv2)%*%matrix(deriv1,ncol=1)
pos.ths<-pre.ths-delta
if (all(abs(delta)<conv)) converged=T
}

deriv2<-makeHessian(ipar,resp,pos.ths,p,sigma,D=D,add.sigma=T)

theta<-as.vector(round(pos.ths,4))
SE<-round(sqrt(abs(diag(solve(deriv2)))),4)

return(list(theta=theta,SE=SE,Hessian=deriv2,niter=iter))
}

