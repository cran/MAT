makeFI <-
function(ipar,th,p,sigma,D=1.7,add.sigma=T) {
#returns a p-by-p Fisher information matrix at given thetas for a group of items

ni<-nrow(ipar)

if (p != length(th)) stop("th and p are non-conforming")


FI<-matrix(0,p,p) #initialize Fisher Information matrix

for (i in 1:ni) {
FI<-FI + calcFI(ipar[i,],th,p,D=D)
}

if (add.sigma) FI<-FI + solve(sigma)
return(FI)
}

