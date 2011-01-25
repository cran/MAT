dLL <-
function(ipar,resp,th,p,sigma,D=1.7,Bayesian=T) {
#returns first-order partial derivative of the log likelihood function
if (p != length(th)) stop("th and p are non-conforming")
ni<-nrow(ipar)

if (p==1) sigma.inv=1
else sigma.inv<-solve(sigma)

dll<-rep(0,p)

aa<-matrix(ipar[,1:p],ncol=p)
dd<-ipar[,p+1]
cc<-ipar[,p+2]

for (i in 1:ni) {
u<-resp[[i]]
if (u==1 || u==0) {
a<-matrix(aa[i,],nrow=1)
d<-dd[i] #d=-a*b
c<-cc[i] 

P<-c + (1-c)/(1+exp(-D*(a%*%th + d)))

dll<-dll + as.vector(a)*(P-c)*(u-P)/((1-c)*P)
}
}

if (Bayesian) {
for (h in 1:p) {
w<-matrix(0,ncol=p)
w[h]<-1
dll[h]<-dll[h] - w%*%sigma.inv%*%th
}
}

dll<-D*dll
return(dll)
}

