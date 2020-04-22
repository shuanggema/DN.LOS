#####################################################Functions used

#Conditional distribution
cond.prob.G<-function(G,H,K,v,y,b){#omage1 is scale, omega2 is vector
  g<-G[b,b]+2*G[b,-b]%*%v[-b]+H[b,-b]%*%y[-b]
  h<-H[b,b]+t(H[-b,b])%*%v[-b]-K[b,-b]%*%y[-b]
  k<-K[b,b]
  prec<-k
  mu<-h/k
  c<-h^2/2/k+g-log(k/2/pi)/2
  p<-1/(1+exp(-c))
  return(ifelse(runif(1)<p,rnorm(1,mu,sqrt(1/prec)),0))
}

#Gibbs sampling
hur.gibb.G<-function(n,G,H,K,thre){
  p.dim<-dim(G)[1]
  simu.data.y<-matrix(0,nrow=thre+n,ncol=p.dim)
  
  y.simu<-rep(0,p.dim)
  for (i in 1:(thre+n)){
    for (l in 1:p.dim){
      y.trans<-cond.prob.G(G,H,K,v=ifelse(y.simu!=0,1,0),y=y.simu,b=l)
      y.simu[l]<-y.trans
    }
    simu.data.y[i,]<- y.simu
  }
  simu.data.y<-simu.data.y[(thre+1):(thre+n),]
  return(simu.data.y)
}

#Selection sampling
hur.sel.G<-function(n,mu,K,a,b){
  library(MASS)
  library(pracma)
  y<-mvrnorm(n,rep(mu,ncol(K)),inv(K))
  v<-y
  for(i in 1:n){
    for(j in 1:ncol(K)){
      x<-a[j]+b*y[i,j]
      p<-1/(1+exp(-x))
      v[i,j]<-rbinom(1,1,p)
    }
  }
  return(v*y)
}

#Generate band matrices for chain graph topology
Omega.ge<-function(p.dim,diag=1,offdiag=0.5,size=1){
  
  omega<-matrix(0,ncol=p.dim,nrow=p.dim)
  
  diag(omega)<-diag
  for(i in 2:p.dim){
    for (j in max((i-size),1):(i-1)){
      omega[i,j]<-omega[j,i]<-offdiag
    }
  }
  
  return(omega)
}

##Generate matrices for random graph topology
Theta.ge<-function(p.dim,diag,offdiag,sparsity){
  
  theta<-matrix(0,p.dim,p.dim)
  
  diag(theta)<-diag
  for (i in 1:(p.dim-1)){
    for(j in (i+1):p.dim){
      if(runif(1,0,1)<sparsity){
        theta[i,j]<-theta[j,i]<-offdiag
      }
    }
  }
  theta
}

#Contaminate non-zero values with a t8-distributed noise scaled to match the standard deviation of each coordinate
contam <- function(x, df=8, tol=0){
  for(i in seq_len(ncol(x))){
    xi = x[,i]
    xpos <- abs(xi)>tol
    sx <- sd(xi[xpos])
    if(is.na(sx)){sx<-0}
    x[xpos,i] <- xi[xpos] + rt(sum(xpos), df=df)*sx
  }
  x
}


##########################################################Simulation

#Generate the interaction matrices
G.real<-Omega.ge(100,1,-2,1)
H.real<-Omega.ge(100,1,0,1)
K.real<-Omega.ge(100,1,0,1)


eK <- min(eigen(K.real, only.values=TRUE)$values-.1, 0)
diag(K.real) <- diag(K.real)-diag(K.real)*eK #Keep it PD 


#Gibbs sampling
a=hur.gibb.G(50000,G.real,H.real,K.real,200000)
x=c(1:50000)
y=a[x%%5==1,]


#t8-distributed noise
y<-contam(y)


#Selection sampling
y<-hur.sel.G(10000,0,K.real,runif(100,-2,0),1)










