if (!'plyr' %in% installed.packages()[,1]) install.packages(plyr)   # check if plyr installed
library(plyr)
 
# setup
mc.seed    = 100000
burn.in    = 2000   
step.size  = 4
n.mcmc     = 3000*step.size + burn.in  

g = function(arg){c(arg[1],arg[1]/(1+arg[1]^4),cos(1.2*arg[2]))}

fullxt = function(t,yt,xtm1,xt,xtp1,theta,sig,tau){
  zm1 = g(c(xtm1,t-1))
  z   = g(c(xt,t))
  full = dnorm(xt,sum(zm1*theta),tau,log=TRUE)+
         dnorm(xtp1,sum(z*theta),tau,log=TRUE)+
         dnorm(yt,xt^2/20,sig,log=TRUE)
  return(full)
}

fullxn = function(t,yt,xtm1,xt,theta,sig,tau){
  zm1  = g(c(xtm1,t-1))
  z    = g(c(xt,t))
  full = dnorm(xt,sum(zm1*theta),tau,log=TRUE)+
         dnorm(yt,xt^2/20,sig,log=TRUE)
  return(full)
}

pr <- progress_text()           # displays progress
pr$init(n.mcmc)

draw.x = function(y,x,x0,theta,sig,tau,v){
  x1  = rnorm(1,x[1],v)
  num = fullxt(1,y[1],x0,x1  ,x[2],theta,sig,tau)
  den = fullxt(1,y[1],x0,x[1],x[2],theta,sig,tau)
  if (log(runif(1))<(num-den)){x[1]=x1}
  xn  = rnorm(1,x[n],v)
  num = fullxn(n,y[n],x[n-1],xn  ,theta,sig,tau)
  den = fullxn(n,y[n],x[n-1],x[n],theta,sig,tau)
  if (log(runif(1))<(num-den)){x[n]=xn} 
  for (t in 2:(n-1)){
    xt  = rnorm(1,x[t],v)
    num = fullxt(t,y[t],x[t-1],xt  ,x[t+1],theta,sig,tau)
    den = fullxt(t,y[t],x[t-1],x[t],x[t+1],theta,sig,tau)
    if (log(runif(1))<(num-den)){x[t]=xt}
  }
  pr$step()
  return(x)
}

fixedpar = function(y,X,b,A,v,lam){
  n     = length(y)
  k     = ncol(X)
  par1  = (v+n)/2
  var   = solve(crossprod(X,X)+A)
  mean  = matrix(var%*%(crossprod(X,y)+crossprod(A,b)),k,1)
  par2  = v*lam + sum((y-crossprod(t(X),mean))^2)
  par2  = (par2 + crossprod(t(crossprod(mean-b,A)),mean-b))/2
  sig2  = 1/rgamma(1,par1,par2)
  var   = var*sig2
  mean  = mean + crossprod(t(chol(var)),rnorm(k))
  return(c(mean,sig2))
}

quant005   = function(x){quantile(x,0.005)}
quant995   = function(x){quantile(x,0.995)}

# Simulating the data
set.seed(mc.seed)
n         = 100
sig2      = 1
tau2      = 10
sig       = sqrt(sig2)
tau       = sqrt(tau2)
theta     = c(0.5,25,8)
ptr       = c(theta,tau2,sig2)
y         = rep(0,n)
x         = rep(0,n)
x0        = 0
Z         = g(c(x0,0))
x[1]      = rnorm(1,sum(Z*theta),tau)
y[1]      = rnorm(1,x[1]^2/20,sig)
for (t in 2:n){
  Z    = g(c(x[t-1],t-1))
  x[t] = rnorm(1,sum(Z*theta),tau)
  y[t] = rnorm(1,x[t]^2/20,sig)
}
xtr = x

# process graphics
par(mfrow=c(2,2), mar=c(2.5,2.5,1,1), mgp=c(1.6,.6,0))
ts.plot(y, xlab="time", ylab=expression(y[~t]))
ts.plot(x, xlab="time", ylab=expression(x[~t]))
plot(x, y, xlab=expression(x[~t]), ylab=expression(y[~t]))
plot.ts(x[-n], x[-1], ylab=expression(x[~t]), xlab=expression(x[~t-1]), cex=.8)


# Prior hyperparameters
m0        = 0.0
C0        = 10.0
n0        = 6
sig20     = 2/3
theta0    = c(0.5,25,8)
V0        = diag(c(0.25,10,4))/tau2
nu0       = 6
tau20     = 20/3
iV0       = diag(1/diag(V0))
n0sig202  = n0*sig20/2
            c((n0*sig20/2)/(n0/2-1),((n0*sig20/2)/(n0/2-1))^2/(n0/2-2),
            (nu0*tau20/2)/(nu0/2-1),((nu0*tau20/2)/(nu0/2-1))^2/(nu0/2-2))

# MCMC setup
v     = 0.5
niter = n.mcmc   
theta = ptr[1:3]
tau2  = ptr[4]
sig2  = ptr[5]
x     = xtr
xs    = NULL
ps    = NULL
n0n2  = (n0+n)/2
for (iter in 1:niter){	
  x     = draw.x(y,x,x0,theta,sig,tau,v)
  sig2  = 1/rgamma(1,n0n2,n0sig202+sum((y-x^2/20)^2)/2)
  sig   = sqrt(sig2)
  W     = cbind(c(x0,x[1:(n-1)]),0:(n-1))
  Z     = t(apply(W,1,g))
  par   = fixedpar(x,Z,theta0,iV0,nu0,tau20)
  theta = par[1:3]
  tau2  = par[4]
  tau   = sqrt(tau2)
  xs    = rbind(xs,x)
  ps    = rbind(ps,c(par,sig2))
}

M0 = burn.in     
M  = n.mcmc-M0  
L  = step.size   

index = seq((M0+1),niter,by=L) 
length(index)

mx = apply(xs[index,],2,median)
lx = apply(xs[index,],2,quant005)
ux = apply(xs[index,],2,quant995)

# display state estimation
dev.new()
par(mar=c(3,3,1.5,1), mgp=c(1.6,.6,0))  
plot(xtr,type="p",pch=20,ylim=range(lx,ux+10,xtr),xlab="time",ylab="",main="State Estimation")
lines(mx,col=1,type="o")
legend("topleft",legend=c("True","Posterior Median","99%CI"),col=c(1,1,gray(.5, alpha=.4)),
       lty=c(0,1,1), pch=c(20,1,-1), lwd=c(1,1,5), cex=.8)
 xx=c(time(xtr), rev(time(xtr)))
 yy=c(lx, rev(ux))
polygon(xx, yy, border=NA, col=gray(.5, alpha = .4))


# display parameter estimation
names = c(expression(alpha), expression(beta), expression(gamma),
          expression(sigma[w]^2), expression(sigma[v]^2))

dev.new()
par(mfrow=c(3,5), mar=c(2.75,2.75,2,1), mgp=c(1.6,.6,0))
for (i in 1:5){
  plot(1:length(index),ps[index,i],type="l",xlab="iteration",ylab="")
  title(names[i], cex=1.5)
  abline(h=ptr[i],col=2,lwd=4)
}
for (i in 1:5)
  acf(ps[index,i],main="")
for (i in 1:5){
  hist(ps[index,i],prob=TRUE, xlab="",ylab="",main="")
  lines(density(ps[index,i],adjust=2))
   lp = apply(as.matrix(ps[index,i]),2,quant005)
   up = apply(as.matrix(ps[index,i]),2,quant995)
  abline(v=c(ptr[i],lp,up),col=c(2,4,4),lwd=c(4,1,1), lty=c(1,2,2))
}
