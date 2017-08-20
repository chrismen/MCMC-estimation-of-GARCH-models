


gc()
rm()

mydata = read.csv("simulated_return.csv")
T1=10000
N=20000
z=rep(0,T1)
y=mydata$y
alpha0=0.1
alpha1=0.4
beta1=0.4

k=1
dd=data.frame(k,alpha0,alpha1,beta1)
write.csv(dd, file = "R-mcmc-time-series.csv",row.names=FALSE)
z[1]=y[1]*y[1];
for ( j in 2:T1)
{
  z[j]=alpha0 + alpha1 * y[j-1]*y[j-1]+ beta1*z[j-1];
}
for (k in 1:N)
{
  left=0.0;
  right=30.0;
  for ( j in 2:T1)
  {
    u1=runif(1, min = 0, max = 1)
    b1=u1/sqrt( alpha0 + alpha1*y[j-1]*y[j-1] +beta1*z[j-1])
    c1_right =1.0/b1/b1 -alpha1*y[j-1]*y[j-1] -beta1*z[j-1]
    u1=runif(1, min = 0, max = 1)
    b1=u1*exp(  -y[j]*y[j]/2.0/( alpha0 + alpha1*y[j-1]*y[j-1] + beta1*z[j-1]))
    c1_left=-y[j]*y[j]/2.0/log(b1) -alpha1*y[j-1]*y[j-1]- beta1*z[j-1]
    left=max(left,c1_left)
    right=min(c1_right,right)
  }
  u=runif(1, min = 0, max = 1)
  alpha0=left + (right-left)*u
  left=0.0;
  right=30.0;
  for ( j in 2:T1)
  {
    u1=runif(1, min = 0, max = 1)
    b1=u1/sqrt( alpha0 + alpha1*y[j-1]*y[j-1]  +beta1*z[j-1]  )
    c1_right =(  1.0/b1/b1 -alpha0 -beta1*z[j-1])/y[j-1]/y[j-1]
    u1=runif(1, min = 0, max = 1)
    b1=u1*exp(  -y[j]*y[j]/2.0/(alpha0+alpha1*y[j-1]*y[j-1]+beta1*z[j-1] )  )
    c1_left=  ( -y[j]*y[j]/2.0/log(b1) -alpha0 -beta1*z[j-1])/y[j-1]/y[j-1]
    left=max(left,c1_left)
    right=min(c1_right,right)
  }
  right=min(1.0-beta1,right)
  u=runif(1, min = 0, max = 1)
  alpha1=left + (right-left)*u;
  left=0.0;
  right=30.0;
  for ( j in 2:T1)
  {
    u1=runif(1, min = 0, max = 1)
    b1=u1/sqrt( alpha0+alpha1*y[j-1]*y[j-1]  +beta1*z[j-1] )
    c1_right =( 1.0/b1/b1 -alpha0 -alpha1*y[j-1]*y[j-1])/z[j-1]
    u1=runif(1, min = 0, max = 1)
    b1=u1*exp(  -y[j]*y[j]/2.0/(alpha0+alpha1*y[j-1]*y[j-1]+beta1*z[j-1] )  )
    c1_left=  (  -y[j]*y[j]/2.0/log(b1) -alpha0 -alpha1*y[j-1]*y[j-1])/z[j-1]
    left=max(left,c1_left)
    right=min(c1_right,right)
  }
  right=min(1.0-alpha1,right)
  u=runif(1, min = 0, max = 1)
  beta1=left + (right-left)*u
  z[1]=y[1]*y[1];
  for ( j in 2:T1)
  {
    z[j]=alpha0 + alpha1 * y[j-1]*y[j-1]+ beta1*z[j-1]
  }
  dd=data.frame(k,alpha0,alpha1,beta1)
  write.table(dd,"R-mcmc-time-series.csv", row.names=F,na="NA",append=T, sep=",", col.names=F);
}


