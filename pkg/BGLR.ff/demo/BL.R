#This example will run a standard Bayesian LASSO

rm(list=ls())
setwd(tempdir())

data(wheat)
set.seed(12345)
varB<-0.5*(1/sum(apply(X=wheat.X,MARGIN=2,FUN=var)))
b0<-rnorm(n=1279,sd=sqrt(varB))
signal<-wheat.X%*%b0
error<-rnorm(599,sd=sqrt(0.5))
y<-100+signal+error
 	
nIter=100;
burnIn=50;
thin=10;
saveAt='';
S0=NULL;
weights=NULL;
R2=0.5;
ETA<-list(list(X=wheat.X,model='BL'))
  
fit_BL=BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
plot(fit_BL$yHat,y)
