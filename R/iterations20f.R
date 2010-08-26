nonnegative<-FALSE
print("loading...")
load("spectra.R")
load("noise.R")
load("mask.R")
load("lamda_tel.R")
N<-nrow(spectra)
M<-length(lamda_tel)

source("functions.R")

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

for (Kminus1 in 1:3) {

for (epsilon in c(100.,30.,10.,3.,1.)) {
#epsilon<-1.

## initialize
print("initializing...")
if (nonnegative){
## make the data non-negative, as weakly as possible
spectra200[spectra200<0]<-0.01*noise200[spectra200<0]
K<-1+Kminus1
tiny<-0.01*mean(abs(spectra200))
g<-as.matrix(abs(spectra200[1:K,])+tiny)
g<-g/normbase(g)
a<-outer(sqrt(rowMeans(spectra200^2)),rep(1./K,K))
epsilon<-0
onomadir<-as.character(paste("K",K,"M",M,"N",N,"_nn",sep=""))
}else{
K<-1+Kminus1 
onomabase<-as.character(paste(Kminus1,"pc_basis.R",sep=""))
load(onomabase)
K<-nrow(pc_basis)
g<-pc_basis/normbase(pc_basis)
g<-g/normbase(g)
onomadir<-as.character(paste("K",K,"M",M,"N",N,"ep",epsilon,sep=""))
}
system(paste("mkdir -p",onomadir,sep=" "))
onomabase0<-as.character(paste(onomadir,"/g_0.R",sep=""))
save(g,file=onomabase0)

## get everything into matrix form for iterations
spectra200<-as.matrix(spectra200)
invvar200<-as.matrix(1.0/noise200^2)

print("iterating...")
if(nonnegative){
  n_iter<-200
}else{
  n_iter<-20
}
xi1wn<-vector(length=n_iter)
xk1wn<-vector(length=n_iter)
for (m in 1:n_iter){
print(m)

if(nonnegative){
a<-astepnn(a,spectra200,invvar200,g)
}else{
a<-astep(spectra200,noise200,g)
}
xi1wn[m]<-badness(a,g,spectra200,noise200,0.0)

if(nonnegative){
g<-gstepnn(g,spectra200,invvar200,a)
}else{
g<-gstep(spectra200,noise200,a,g,epsilon)
}
xk1wn[m]<-badness(a,g,spectra200,noise200,epsilon)

norm<-normbase(g)
g<-g/norm
a<-t(t(a)*norm)
print(c(xi1wn[m],xk1wn[m],range(a),range(g)))

## Save coefficients from the 1st fitting (a), the new base from the 2nd f (base), the sum of the x2 of all the spectra (xi, xi3 - 1st f), the sum of the the x2 for all the lamda (xk, xk3 - 2nd f) and the x2 for every spectrum (xi5 -1st f)
epssuffix<-as.character(paste("_",log10(epsilon),".R",sep=""))
suffix<-as.character(paste("_",m,epssuffix,sep=""))
onoma<-as.character(paste(onomadir,"/a",suffix,sep=""))
save(a,file=onoma)
onoma1<-as.character(paste(onomadir,"/g",suffix,sep=""))
save(g,file=onoma1)
}

onoma11<-as.character(paste(onomadir,"/xi1wn",epssuffix,sep=""))
save(xi1wn,file=onoma11)
onoma12<-as.character(paste(onomadir,"/xk1wn",epssuffix,sep=""))
save(xk1wn,file=onoma12)
}
}
