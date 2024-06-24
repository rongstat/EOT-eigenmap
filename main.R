#######################################################
###     Entropic Optimal Transport (EOT) Eigenmap
#######################################################

library(Rfast)
library(RSpectra)

# X: first data matrix (pxn1) to be aligned.
# Y: second data matrix (pxn2) to be aligned.
# bdw: numeric (default NULL); user-specified bandwith parameter for the kernel matrix. If set as NULL, a data-driven procedure described in the paper will be used to determine the bandwidth.
# t: integer (default 0); the diffucion steps.
# q: integer (default NULL); the output embedding dimension. If set as NULL, a data-driven procedure described in the paper will be used to determine its value.
# tau: numeric (default 0.02); convergence threshold for the Sinkhorn-Knopp algorithm.
# output: a list containing the EOT low-dimensional embedded dataset.

EOTmap <- function(X, Y, bdw = NULL, t=0, q=NULL, tau=0.02){
  
  
  dist.mat = Dist(cbind(X, Y))
  cross.dist = as.matrix(dist.mat)[1:dim(X)[2], (dim(X)[2]+1):(dim(X)[2]+dim(Y)[2])]
  rm(dist.mat)
  if(is.null(bdw)){
    K.mat = exp(-cross.dist^2/quantile(cross.dist,0.5)^2) 
  }else{
    K.mat = exp(-cross.dist^2/bdw^2) 
  }
  ##Sinkhorn-Knopp algorithm
  u=rep(1,dim(K.mat)[1])
  v=rep(1,dim(K.mat)[2])
  r.s = rep(dim(K.mat)[2],dim(K.mat)[1])
  c.s = rep(dim(K.mat)[1],dim(K.mat)[2])
  while(max(abs(diag(u) %*% K.mat %*% v-r.s))>tau | max(abs( t(u)%*% K.mat %*% diag(v)-c.s))>tau){
    v <- c(c.s/(t(u) %*% K.mat))
    u <- c(r.s/(K.mat %*% v))
  }
  print("EOT plan obtained!")
  K.mat.OT = diag(u) %*% K.mat %*% diag(v)
  if(is.null(q)){
    prop.svd = svds(K.mat.OT, k=80)
    r=max(which(prop.svd$d[1:79]/prop.svd$d[2:80]>1.02))
  }else{
    prop.svd = svds(K.mat.OT, k=q+1)
    r=q
  }
  print("EOT embedding obtained!")
  out.emb1 = sqrt(dim(X)[2])*prop.svd$u[,1:r+1] %*% diag(prop.svd$d[2:(r+1)]^t)
  out.emb2 = sqrt(dim(Y)[2])*prop.svd$v[,1:r+1] %*% diag(prop.svd$d[2:(r+1)]^t)
  
  return(list(X.embed=out.emb1, Y.embed=out.emb2))
}
