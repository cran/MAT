SCORE <-
function(ipar,resp,p,sigma,maxIter=30,conv=0.001,D=1.7,Fisher=T) {
#calculates MAP theta and SE estimates for a group
ni<-nrow(ipar) #number of items
nj<-nrow(resp) #number of persons
TH<-matrix(nrow=nj,ncol=p)
SE<-matrix(nrow=nj,ncol=p)

for (j in 1:nj) {
pos.ths<-numeric(p) 
iter<-0
converged=F
rs<-resp[j,]

while (iter<maxIter && !converged) {
iter<-iter+1
pre.ths<-pos.ths

deriv1<-dLL(ipar,resp[j,],pre.ths,p,sigma,D=D,Bayesian=T)

if (Fisher) deriv2<- -makeFI(ipar,pre.ths,p,sigma,D=D,add.sigma=T)
else deriv2<-makeHessian(ipar,resp[j,],pre.ths,p,sigma,D=D,add.sigma=T)

delta<-solve(deriv2)%*%matrix(deriv1,ncol=1)
pos.ths<-pre.ths-delta
if (all(abs(delta)<conv)) converged=T
}

deriv2<-makeHessian(ipar,resp[j,],pos.ths,p,sigma,D=D,add.sigma=T)

TH[j,]<-pos.ths
SE[j,]<-sqrt(abs(diag(solve(deriv2))))
}
return(list(theta=round(TH,4),SE=round(SE,4)))
}

