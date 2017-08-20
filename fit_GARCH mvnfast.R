

library("microbenchmark")
library("mvtnorm")
library("mvnfast")
library("RhpcBLASctl")
blas_set_num_threads(1)

mydata = read.csv("R-mcmc-time-series.csv", header= TRUE)
aa=data.matrix(mydata[10001:20000,2:4])
mcov1=cov(aa)*8

T=10000
z=rep(0,T)
alpha0=0.1
alpha1=0.2
beta1=0.14
mu=c(alpha0, alpha1, beta1)
mu1=rmvn(1, mu, mcov1, ncores = 2)
mydata = read.csv("simulated_return.csv")
y=mydata$y
dd=data.frame(1,alpha0,alpha1,beta1)
write.table( dd, file = "mcmc_mnormal fast.csv",sep=",",row.names=FALSE, col.names=FALSE)
z[1]=y[1]*y[1];
for ( j in 2:T)
{
  z[j]=alpha0 + alpha1 * y[j-1]*y[j-1]+ beta1*z[j-1];
}
N=20000
tau=0.1
mcov1=tau*mcov1
for (k in 1:N)
{
  mu1= rmvn(1, mu, mcov1, ncores = 2)
  likelihood=0
  talpha0=mu1[1]
  talpha1=mu1[2]
  tbeta1=mu1[3]
  for ( j in 2:T)
  {
    likelihood=likelihood -y[j]*y[j]/2.0/(talpha0+talpha1*y[j-1]*y[j-1]+tbeta1*z[j-1] ) +y[j]*y[j]/2.0/(alpha0+alpha1*y[j-1]*y[j-1]+beta1*z[j-1] ) 
    likelihood=likelihood +log(sqrt(  alpha0+alpha1*y[j-1]*y[j-1]  +beta1*z[j-1] ))-log(sqrt(talpha0+talpha1*y[j-1]*y[j-1]  +tbeta1*z[j-1] ))
  }
  u=runif(1, min=0, max=1)
  alpha=min(1, exp(likelihood))
  if (u<alpha){
    alpha0=mu1[1]
    alpha1=mu1[2]
    beta1=mu1[3]
    mu=mu1
  }
  for ( j in 2:T)
  {
    z[j]=alpha0 + alpha1 * y[j-1]*y[j-1]+ beta1*z[j-1]
  }
  dd=data.frame(k,alpha0,alpha1,beta1)
  write.table(dd,"mcmc_mnormal fast.csv", row.names=F,na="NA",append=T, sep=",", col.names=F);
  print(k)
}




