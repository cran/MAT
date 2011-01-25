calcFI <-
function(ipar,th,p,D=1.7) {
#returns a p-by-p Fisher information matrix at given thetas for a single item

if (p != length(th)) stop("th and p are non-conforming")

a<-matrix(ipar[1:p],nrow=1)
d<-ipar[p+1] #d=-a*b
c<-ipar[p+2]

P<-c + (1-c)/(1+exp(-D*(a%*%th + d))) 

cf<-as.numeric((1-P)*(P-c)^2/(P*(1-c)^2))

FI<-matrix(a,ncol=1)%*%matrix(a,nrow=1)

FI<-D^2*cf*FI

return(FI)
}

