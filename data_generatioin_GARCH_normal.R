

gc()
rm()

x=rep(0,2000)
y=rep(0,2000)
alpha0=0.1
alpha1=0.4
beta1=0.4
x1=sqrt( alpha0/(1-alpha1-beta1))
y[1]=x1*rnorm(1, mean=0,sd=1)
k1=20000
T1=10000
for(i in 1:k1)
{
  x2=alpha0 + alpha1*y[1]^2 + beta1*x1
  y[1]=sqrt(x2)*rnorm(1, mean=0,sd=1)
  x1=x2
  print(y[1])
}
x2=alpha0 + alpha1*y[1]^2 + beta1*x1
y[1]=sqrt(x2)*rnorm(1, mean=0,sd=1)
x1=x2
x[1]=x2
for (j in 2:T1)
{
  x2=alpha0 + alpha1*y[j-1]^2+ beta1*x1
  y[j]=sqrt(x2)*rnorm(1, mean=0, sd=1)
  x1=x2;
  x[j]=x2;
  print(y[j])
}
mydata=data.frame(y,x)
colnames(mydata) <- c("y", "x")
write.csv(mydata, file = "simulated_return.csv",row.names=FALSE)

dy <- data.frame(y)    
write.table(dy, "C-simulated-returns.txt",row.names = FALSE, col.names = FALSE)







