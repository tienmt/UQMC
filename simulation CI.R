P_Omega = function(a,entri){
  a[entri] = 0
  return(a)
}
P_Omega <- compiler::cmpfun(P_Omega)
vij.get = function(i,j,alpha,Xd,Yd...){
  v.ij = sig*( Xd[i,] %*%solve(crossprod(Xd))%*% Xd[i,] + Yd[j,] %*%solve(crossprod(Yd))%*% Yd[j,])
  v.ij = c(x[i,j] - qnorm(1-alpha/2)*sqrt(v.ij),x[i,j] + qnorm(1-alpha/2)*sqrt(v.ij))
  return(v.ij)
}
vij.get <- compiler::cmpfun(vij.get)
alpha = 0.05
library(softImpute)
library(RSpectra)
p = 1000
n = 100

J = 5
np = n*p
missfrac = 0.8
sig = 1
lamda = 2.5*sig*sqrt(n*p)
rk = J
x = matrix(rnorm(n*J),n,J)%*%matrix(rnorm(J*p),J,p) #+ 0.1*matrix(rnorm(n*50),nr=n)%*%matrix(rnorm(50*p),nc=p) 

Xbay = Xfbay = ci= x11= list()
als.mse = als.nmse = als.pre =
  db.mse = db.nmse = db.pre = 
  fbay.mse = fbay.nmse = fbay.pre = 
  bay.mse = bay.nmse = bay.pre =  c()
for(s in 1:50){
  imiss = sample(np,np*missfrac,replace=FALSE)
  xna = x + matrix(rnorm(np),nr = n,nc = p)
  xna[imiss] = NA
  x11[[s]] = x[1,1]
  # ALS method
  fit1 = softImpute(xna, maxit = 1000,type = 'als',rank.max = J)
  ximp = complete(xna,fit1)
  # de-biased method
  tam = svds(ximp - P_Omega(ximp - xna,imiss),rk)
  lamda = 2.5*sig*sqrt(n*p)
  Xd = tam$u[,1:rk] %*% diag( sqrt(tam$d[1:rk] + lamda) )
  Yd = tam$v[,1:rk]%*% diag( sqrt(tam$d[1:rk] + lamda) )
  ci[[s]] = vij.get(1,1,alpha,Xd,Yd)
  x.db = tam$u[,1:rk] %*% diag(tam$d[1:rk]) %*% t(tam$v[,1:rk])
  db.mse[s] = mean((x.db - x)^2)
  
  #Gibbs for fix k
  I1 = row(xna)[!is.na(xna)]
  I2 = col(xna)[!is.na(xna)]
  Y = xna[!is.na(xna)]
  nsample =  (1-missfrac)*n*p
  obs = Y
  sigma = 1
  
  k = J#min(n,p)/2 #rank of U,V
  lambda = 1/(4*sigma^2)
  Mstep = matrix(1 ,nr=n,nc=k)
  Nstep = matrix(1,nr=p,nc=k)
  gamma = rep(1,k)
  Xmean = matrix(data=0,nr=n,nc=p)
  L2 = lambda
  Nmcmc = 500
  burnin = 100
  Xfbay[[s]] = list()
  datalength = as.vector(1:nsample)
  for (step in 1:(Nmcmc+burnin)){  
    # update M[i,j]
    for (i in 1:n) {
      seti = datalength[I1==i]
      seti = seti[!is.na(seti)]
      for (j in 1:k){
        Msteptrouj = Mstep[i,]
        Msteptrouj[j] = 0
        V = (1/gamma[j]) + L2*sum(Nstep[I2[seti],j]^2)
        D = sum(L2*(obs[seti] - Msteptrouj%*%t(Nstep[I2[seti],]) )*Nstep[I2[seti],j])
        Mstep[i,j] = rnorm(1,D/sqrt(V),1)/sqrt(V)
      }
    }
    # update N[i,j]
    for (i in 1:p) {
      seti = (1:nsample)[I2==i]
      for (j in 1:k){
        Nsteptrouj = Nstep[i,]
        Nsteptrouj[j] = 0
        V = (1/gamma[j]) + L2*sum(Mstep[I1[seti],j]^2)
        D = sum(L2*(obs[seti] - Nsteptrouj%*%t(Mstep[I1[seti],]))*Mstep[I1[seti],j])
        Nstep[i,j] = rnorm(1,D/sqrt(V),1)/sqrt(V)
      }
    }
    l1 <- list(Xmean*(1-1/(step-burnin)), Mstep%*%t(Nstep)/(step-burnin) )
    res <- Reduce(`+`, lapply(l1, function(x)replace(x, is.na(x), 0)))
    Xmean = res*NA^!Reduce(`+`, lapply(l1, function(x) !is.na(x)))
    if(step> burnin) Xfbay[[s]][[step]] =  Mstep%*%t(Nstep)
  }
  fbay.mse[s] = mean((Xmean - x )^2)
  
  # full Bayes Gibbs sampler
  k = 10#min(n,p)/2 #rank of U,V
  a = 1
  b = 1/100
  lambda = 1/(4*sigma^2)
  Mstep = matrix(1 ,nr=n,nc=k)
  Nstep = matrix(1,nr=p,nc=k)
  gamma = rep(b/a,k)
  Xm = matrix(data=0,nr=n,nc=p)
  L2 = 2*lambda
  Xbay[[s]] = list()
  Nmcmc = 500
  burnin = 100
  datalength = as.vector(1:nsample)
  for (step in 1:(Nmcmc+burnin)){  
    # update M[i,j]
    for (i in 1:n) {
      seti = datalength[I1==i]
      seti = seti[!is.na(seti)]
      for (j in 1:k){
        Msteptrouj = Mstep[i,]
        Msteptrouj[j] = 0
        V = (1/gamma[j]) + L2*sum(Nstep[I2[seti],j]^2)
        D = sum(L2*(obs[seti] - Msteptrouj%*%t(Nstep[I2[seti],]) )*Nstep[I2[seti],j])
        Mstep[i,j] = rnorm(1,D/sqrt(V),1)/sqrt(V)
      }
    }
    # update N[i,j]
    for (i in 1:p) {
      seti = (1:nsample)[I2==i]
      for (j in 1:k){
        Nsteptrouj = Nstep[i,]
        Nsteptrouj[j] = 0
        V = (1/gamma[j]) + L2*sum(Mstep[I1[seti],j]^2)
        D = sum(L2*(obs[seti] - Nsteptrouj%*%t(Mstep[I1[seti],]))*Mstep[I1[seti],j])
        Nstep[i,j] = rnorm(1,D/sqrt(V),1)/sqrt(V)
      }
    }
    # update gamma
    for (j in 1:k) gamma[j] = 1/rgamma(1,a+(n+p)/2,b+(sum(Mstep[,j]^2)+sum(Nstep[,j]^2))/2)
    l1 <- list(Xm*(1-1/(step-burnin)), Mstep%*%t(Nstep)/(step-burnin) )
    res <- Reduce(`+`, lapply(l1, function(x)replace(x, is.na(x), 0)))
    Xm = res*NA^!Reduce(`+`, lapply(l1, function(x) !is.na(x)))
    if(step> burnin) Xbay[[s]][[step]]  = Mstep%*%t(Nstep)
  }
  bay.mse[s] = mean((Xm - x )^2)
  print(s)
}

x11.Xbay = unlist(lapply(Xbay, lapply,function(x) x[1,1] ))
x11.Xfbay = unlist(lapply(Xfbay, lapply,function(x) x[1,1] ))
save(x11.Xfbay,ci,x11,db.mse,fbay.mse,bay.mse,x11.Xbay,
     file = '/data2/thetm/matrixCompletion/UQ.r5p1000mis02.rda')

