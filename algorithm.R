###############################Line 2-27 runs before estimating each variable seperately

G.initial=1;H.initial=1;K.initial=1;thr=0.00001;t1=0.01;t2=0.01;t3=0.01;lambda=0.01
n<-dim(y)[1]
p<-dim(y)[2]
node<-1:p
G<-matrix(0,nrow=p,ncol=p)
diag(G)<-G.initial
H<-matrix(0,nrow=p,ncol=p)
diag(H)<-H.initial
K<-matrix(0,nrow=p,ncol=p)
diag(K)<-K.initial

G.final<-G
H1.final<-H
H2.final<-H
K.final<-K
v<-ifelse(y!=0,1,0)
v<-as.matrix(v)
y<-as.matrix(y)

v.v<-t(v)%*%v
v.y<-t(v)%*%y
y.v<-t(v.y)
y.y<-t(y)%*%y
y.sum<-apply(y,2,sum)
y2.sum<-apply(y^2,2,sum)
v.sum<-apply(v,2,sum)

############################Line 37-135 is to estimate the parameters of variable iterf-itert

library(foreach)
library(doParallel)

numCores <- detectCores()
registerDoParallel(numCores)


iterf<-1
itert<-p


q1<-Sys.time()

para.list <- foreach (j = iterf:itert) %dopar% {
g.grad<-rep(0,p)
h1.grad<-rep(0,p)
h2.grad<-rep(0,p)
k.grad<-rep(0,p)
g.j<-G[,j]
h1.j<-H[j,]
h2.j<-H[,j]
k.j<-K[,j]


l=0


while(1){
  l<-l+1
  para.original<-c(g.j,h1.j,h2.j,k.j)
  
  gj<-g.j[j]+2*v[,-j]%*%g.j[-j]+y[,-j]%*%h1.j[-j]
  hj<-h1.j[j]+v[,-j]%*%h2.j[-j]-y[,-j]%*%k.j[-j]
  kj<-rep(k.j[j],n)
  wj<-gj+hj^2/kj/2-log(kj/2/pi)/2
  Nj<-1/(1+exp(-wj))
  Q<-Nj*hj/kj
  P<-(hj^2/2/kj^2+1/2/kj)*Nj
  
  g.grad[j]<-v.sum[j]/n-sum(Nj)/n
  g.grad[-j]<-2*v.v[-j,j]/n-2*t(v[,-j])%*%Nj/n
  h1.grad[j]<-y.sum[j]/n-sum(Q)/n
  h1.grad[-j]<-y.v[-j,j]/n-t(y[,-j])%*%Nj/n
  h2.grad[-j]<-v.y[-j,j]/n-t(v[,-j])%*%Q/n
  k.grad[j]<--0.5*y2.sum[j]/n+sum(P)/n
  k.grad[-j]<--y.y[-j,j]/n+t(y[,-j])%*%Q/n
  
  g.j[j]<-g.j[j]+t1*g.grad[j]
  h1.j[j]<-h2.j[j]<-h1.j[j]+t2*h1.grad[j]
  k.j[j]<-k.j[j]+t3*k.grad[j]
  for (k in node[-j]){
    g.trans<-g.j[k]+t1*g.grad[k]
    if(g.trans<(-lambda*t1)){
      g.j[k]<-g.trans+lambda*t1
    }else if(g.trans>lambda*t1){
      g.j[k]<-g.trans-lambda*t1
    }else{g.j[k]<-0}
    
    h1.trans<-h1.j[k]+t2*h1.grad[k]
    if(h1.trans<(-lambda*t2)){
      h1.j[k]<-h1.trans+lambda*t2
    }else if(h1.trans>lambda*t2){
      h1.j[k]<-h1.trans-lambda*t2
    }else{h1.j[k]<-0}
    
    h2.trans<-h2.j[k]+t2*h2.grad[k]
    if(h2.trans<(-lambda*t2)){
      h2.j[k]<-h2.trans+lambda*t2
    }else if(h2.trans>lambda*t2){
      h2.j[k]<-h2.trans-lambda*t2
    }else{h2.j[k]<-0}
    
    k.trans<-k.j[k]+t3*k.grad[k]
    if(k.trans<(-lambda*t3)){
      k.j[k]<-k.trans+lambda*t3
    }else if(k.trans>lambda*t3){
      k.j[k]<-k.trans-lambda*t3
    }else{k.j[k]<-0}
  }
  
  para.update<-c(g.j,h1.j,h2.j,k.j)
  diff<-para.update-para.original
  if(all(abs(diff)<thr)){
    break
  }else{
    next}
}

list(g.j,h1.j,h2.j,k.j)
}


q2<-Sys.time()
q2-q1

stopImplicitCluster()

#######################################Record the estimating result
for (m in iterf:itert){
  G.final[,m]<-para.list[[m-(iterf-1)]][[1]]
  H1.final[m,]<-t(para.list[[m-(iterf-1)]][[2]])
  H2.final[,m]<-para.list[[m-(iterf-1)]][[3]]
  K.final[,m]<-para.list[[m-(iterf-1)]][[4]]
}



