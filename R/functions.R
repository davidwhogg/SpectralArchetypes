## Copyright 2010 by Paraskevi Tsalmantza & David W. Hogg.
## All rights reserved.

## Comments:
## - Synchronize all variable names with documents in
##   http://trac.astrometry.net/browser/trunk/documents/papers/archetypes
## - Keep matrices in same transpose as those documents too:
##   - a is [N,K]
##   - g is [K,M]
##   - data (spectro200) is [N,M]
##   - innvar is either a scalar or else [N,M]
##   - epsilon is a scalar
## - Use of "noise" (rather than "invvar") is deprecated.
## - astep() and gstep() are deprecated.

## Bugs:
## - astep() and gstep() are slow-o-rama.

## compute model
model<- function(a,g) {
return(a%*%g)
}

## compute residuals
resid<- function(a,g,spectra) {
return(spectra-model(a,g))
}

## compute chi (scaled residual)
chi<- function(a,g,spectra,invvar) {
return(resid(a,g,spectra)*sqrt(invvar))
}

## compute penalty for non-smoothness
penalty<- function(g,epsilon) {
M<-ncol(g)
return(epsilon*sum((g[,-M]-g[,2:M])^2))
}

## compute scalar
badness<- function(a,g,spectra,invvar,epsilon){
return(sum(chi(a,g,spectra,invvar)^2)+penalty(g,epsilon))
}

## apply standard component normalization
normbase<- function(g) {
return(sqrt(rowMeans(as.matrix((g)^2))))
}

## find coefficients at fixed component spectra (A step) 
astep<- function(spectra,noise,g) {
N=nrow(spectra)
K=nrow(g)
a<-matrix(0,N,K)
for (i in 1:N) {
c<-as.vector(1/noise[i,]^2)
ex<-lm(spectra[i,]~t(g)+0,weights=c)
a[i,]<-t(ex$coefficients)
}
return(a)
}

## find component spectra at fixed coefficients (G step)
gstep<- function(spectra,noise,a,oldg,epsilon) {
M<-ncol(spectra)
K<-ncol(a)
g<-matrix(0,K,M)
morea<-rbind(diag(rep(1,K)),diag(rep(1,K)))
modelmatrix<-as.matrix(rbind(a,morea))
for (j in 1:M) {
if (j>1){
  epsjm1<-epsilon
  basejm1<-oldg[,j-1]
} else {
  epsjm1<-0.
  basejm1<-rep(0.,K)
}
if (j<M){
  epsjp1<-epsilon
  basejp1<-oldg[,j+1]
} else {
  epsjp1<-0.
  basejp1 <- rep(0.,K)
}
c1<-c(as.vector(1/noise[,j]^2),rep(epsjm1,K),rep(epsjp1,K))
datamatrix<-as.matrix(c(spectra[,j],basejm1,basejp1))
ex1<-lm(datamatrix~modelmatrix+0,weights=c1)
g[,j]<-t(ex1$coefficients)
}
return(g)
}

## gradient-descent update for coefficients fixed component spectra
## following Shewchuk
## http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
astepgd<- function(a,spectra,invvar,g) {
gradient<-((spectra-a%*%g)*invvar)%*%t(g)
Agradient<-((gradient%*%g)*invvar)%*%t(g)
alpha<-sum(gradient*gradient)/sum(gradient*Agradient)
return(a+alpha*gradient)
}

## left-multiply by the horrifying Q-matrix
## (custom sparse matrix operation)
## needed for gstepgd(), and probably other functions
leftQ<- function(g){
M<-ncol(g)
K<-nrow(g)
diff1<-cbind(matrix(0,K,1),g[,1:M-1]-g[,2:M])
diff2<-cbind(g[,2:M]-g[,1:M-1],matrix(0,K,1))
return(diff1+diff2)
}

## gradient-descent update for component spectra at fixed coefficients
## following Shewchuk
## http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
gstepgd<- function(g,spectra,invvar,a,epsilon) {
gradient<-t(a)%*%((spectra-a%*%g)*invvar)+epsilon*leftQ(g)
Agradient<-t(a)%*%((a%*%gradient)*invvar)-epsilon*leftQ(gradient)
alpha<-sum(gradient*gradient)/sum(gradient*Agradient)
return(g+alpha*gradient)
}

## non-negative update for coefficients fixed component spectra
## following Blanton & Roweis
## http://adsabs.harvard.edu/abs/2007AJ....133..734B
astepnn<- function(a,spectra,invvar,g) {
numerator<-(spectra*invvar)%*%t(g)
denominator<-((a%*%g)*invvar)%*%t(g)
a<-a*(numerator/denominator)
return(a)
}

## non-negative update for component spectra at fixed coefficients
## following Blanton & Roweis
## http://adsabs.harvard.edu/abs/2007AJ....133..734B
gstepnn<- function(g,spectra,invvar,a) {
numerator<-t(a)%*%(spectra*invvar)
denominator<-t(a)%*%((a%*%g)*invvar)
g<-g*(numerator/denominator)
return(g)
}
