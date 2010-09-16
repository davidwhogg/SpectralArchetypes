## Copyright 2010 Paraskevi Tsalmantza & David W. Hogg.
## All rights reserved.

nonnegative<-FALSE
print("loading...")
load("spectra.R")
load("noise.R")
load("mask.R")
load("lamda_tel.R")
N<-nrow(spectra)
M<-length(lamda_tel)

source("/home/vivitsal/DAVID/method/galaxies/functions.R")

## Check which spectra have the whole spectrum at the beggining and at the end
all_blue<-which(mask[,1]==0)
all_red<-which(mask[,ncol(mask)]==0)
all_br<-union(all_blue,all_red)

set.seed(8)
index<-sample(setdiff(c(1:N),all_br),(N-M))

spectra1<-spectra[-index,]
spectrat<-spectra[index,]
spectra<-spectra1
noise1<-noise[-index,]
noiset<-noise[index,]
noise<-noise1
noise<-abs(noise) # should be a no-op
mask1<-mask[-index,]
maskt<-mask[index,]
mask<-mask1

all_blue<-which(mask[,1]==0)
all_red<-which(mask[,ncol(mask)]==0)
all_br<-union(all_blue,all_red)

## Keep the the above spectra and the 200 first ones (because this is a test)
#M<-100
spectra200<-as.data.frame(spectra[union(all_br,1:1000),1:M])
noise200<-as.data.frame(noise[union(all_br,1:1000),1:M])
N<-nrow(spectra200)

if (nonnegative){
print("making the data non-negative...")
## make the data non-negative, as weakly as possible
negativedata<-(spectra200<0)
spectra200[negativedata]<-0
}

## get everything into matrix form for iterations
spectra200<-as.matrix(spectra200)
invvar200<-as.matrix(1.0/noise200^2)

Kminus1list<-c(1:15)
epsilonlist<-c(30.,10.,3.,1.)
nK<-length(Kminus1list)
nepsilon<-length(epsilonlist)
xit<-matrix(0.,nepsilon,nK)
Kindex<-0
for (Kminus1 in Kminus1list) {
Kindex<-Kindex+1
K<-Kminus1+1

epsilonindex<-0
for (epsilon in epsilonlist) {
epsilonindex<-epsilonindex+1

## initialize
print("initializing...")
#initialize with Kmeans
g<-kmeans(spectra200,K)$center
g<-g/normbase(g)
a<-outer(sqrt(rowMeans(spectra200^2)),rep(1./K,K))
if (nonnegative){
  for (aiter in 1:128){
    a<-astepnn(a,spectra200,invvar200,g)
  }
}
if (nonnegative){
if (min(spectra200)<0){ stop("spectra200 has negative elements!") }
if (min(a)<0){ stop("a has negative elements!") }
if (min(g)<0){ stop("g has negative elements!") }
epsilon<-0
onomadir<-as.character(paste("K",K,"M",M,"N",N,"_nn",sep=""))
}else{
onomadir<-as.character(paste("K",K,"M",M,"N",N,sep=""))
}
system(paste("mkdir -p",onomadir,sep=" "))
onomabase0<-as.character(paste(onomadir,"/g_0.R",sep=""))
save(g,file=onomabase0)

print("iterating...")
if(nonnegative){
  n_iter<-2048
}else{
  n_iter<-16
}
xi1wn<-vector(length=n_iter)
xk1wn<-vector(length=n_iter)
for (m in 1:n_iter){
print(m)

oldg<-g

if(nonnegative){
a<-astepnn(a,spectra200,invvar200,g)
}else{
a<-astep(spectra200,invvar200,g)
}
xi1wn[m]<-badness(a,g,spectra200,invvar200,0.0)

if(nonnegative){
g<-gstepnn(g,spectra200,invvar200,a)
}else{
g<-gstep(spectra200,invvar200,a,g,epsilon)
}
xk1wn[m]<-badness(a,g,spectra200,invvar200,epsilon)

## reorder and normalize
if(!nonnegative){
  foo<-reorder(a,g)
  a<-foo$a
  g<-foo$g
}
norm<-normbase(g)
g<-g/norm
a<-t(t(a)*norm)
print(c(xi1wn[m],xk1wn[m],range(a),range(g)))

## Save coefficients from the 1st fitting (a), the new base from the 2nd f (base), the sum of the x2 of all the spectra (xi, xi3 - 1st f), the sum of the the x2 for all the lamda (xk, xk3 - 2nd f) and the x2 for every spectrum (xi5 -1st f)
if (m == 2^round(log2(m))){
epssuffix<-as.character(paste("_",log10(epsilon),".R",sep=""))
suffix<-as.character(paste("_",m,epssuffix,sep=""))
onoma<-as.character(paste(onomadir,"/a",suffix,sep=""))
save(a,file=onoma)
onoma1<-as.character(paste(onomadir,"/g",suffix,sep=""))
save(g,file=onoma1)
}
}

## How well does this converged model do on new data?
spectrat<-as.matrix(spectrat[,1:M])
invvart<-as.matrix(1.0/noiset[,1:M]^2)
at<-astep(spectrat,invvart,g)
xit[epsilonindex,Kindex]<-badness(at,g,spectrat,invvart,0.)

onoma11<-as.character(paste(onomadir,"/xi1wn",epssuffix,sep=""))
save(xi1wn,file=onoma11)
onoma12<-as.character(paste(onomadir,"/xk1wn",epssuffix,sep=""))
save(xk1wn,file=onoma12)
}
}
if (nonnegative){
onoma13<-as.character(paste("xit_",(min(Kminus1list)+1),"_",(max(Kminus1list)+1),"_nn.R",sep=""))
  save(xit,file=onoma13)
}else{
onoma14<-as.character(paste("xit_",(min(Kminus1list)+1),"_",(max(Kminus1list)+1),".R",sep=""))
  save(xit,file=onoma14)
}

