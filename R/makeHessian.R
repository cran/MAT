makeHessian <-
function(ipar,resp,th,p,sigma,D=1.7,add.sigma=T) {
#returns a p-by-p Hessian matrix at given thetas for a group of items

ni<-nrow(ipar)

if (p != length(th)) stop("th and p are non-conforming")

H<-matrix(0,p,p) #initialize Hessian matrix

for (i in 1:ni) {
u<-resp[[i]]
if (u==1 || u==0) {
H<-H + calcHessian(ipar[i,],u,th,p,D=D)
}
}

if (add.sigma) H<-H - solve(sigma)
return(H)
}

