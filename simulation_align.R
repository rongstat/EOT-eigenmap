######## simulation for EOT-embedding

library(Rfast)
library(RSpectra)
library(BiocNeighbors)
library(pcaPP)
library(clusterSim)
library(scatterplot3d)
library(ggplot2)

###load data: torus

set.seed(23)
n=5000
u = runif(n, 0, 2*pi)
v = runif(n, 0, 2*pi)
data = cbind( (2+0.8*cos(u))*cos(v),
              (2+0.8*cos(u))*sin(v),
              0.8*sin(u))
data1=data
colfunc <- colorRampPalette(c("blue","forestgreen", "yellow"))
scatterplot3d(x=data[,1], y=data[,2], z=data[,3], angle=35, color = colfunc(n)[rank(-data[,1])],pch=20,
              xlab=NA, ylab=NA, zlab=NA)

r=3

### knn dist

knn_dist <- function(X,Y){
  #knn graphs
  gx = findKmknn(X, k = 50)$index
  gy = findKmknn(Y, k = 50)$index
  
  out=c()
  for(i in 1:dim(gx)[1]){
    out[i]=length(intersect(gx[i,],gy[i,]))/length(union(gx[i,],gy[i,]))
  }
  return(mean(out))
}



### Setting 1


sig.v=(seq(1,8, length.out= 5)) 
time=matrix(ncol=9,nrow=500)
error12=matrix(ncol=9,nrow=500)
error12.out=matrix(ncol=9,nrow=length(sig.v))
error2=matrix(ncol=9,nrow=500)
error2.out=matrix(ncol=9,nrow=length(sig.v))
error1=matrix(ncol=9,nrow=500)
error1.out=matrix(ncol=9,nrow=length(sig.v))
n=600
theta=n^(2/5)

for(pp in 1:length(sig.v)){
  #p=floor(pp.v[pp])
  p=1000
  sig1 = theta*0.05
  
  for(NN in 1:5){
    
    # #setting 1: basic
    data.X1=data[sample(dim(data)[1],n),]
    data.X2=data[sample(dim(data)[1],n),]
    data.Z=matrix(rnorm(n*p,0,sig1), ncol=p)
    data.Y1 = data.Z
    data.Y1[,1:r] = data.Y1[,1:r]+as.matrix(data.X1*theta)
    data.Y1[,1] = data.Y1[,1] + theta*sig.v[pp]
    data.Y1 = data.Y1*3
    
    data.Z=matrix(rnorm(n*p,0,sig1), ncol=p)
    data.Y2 = data.Z
    data.Y2[,1:r] = data.Y2[,1:r]+as.matrix(data.X2*theta)
    
    #various kernel matrices
    data.Y1=t(data.Y1)
    data.Y2=t(data.Y2)
    dist.mat = Dist(t(cbind(data.Y1, data.Y2)))
    cross.dist = as.matrix(dist.mat)[1:dim(data.Y1)[2], (dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]
    
    K.mat = exp(-cross.dist^2/quantile(cross.dist,0.5)^2) 
    K.mat.rl =  K.mat
    
    ##Sinkhorn-Knopp algorithm
    u=rep(1,dim(K.mat)[1])
    v=rep(1,dim(K.mat)[2])
    r.s = rep(dim(K.mat)[2],dim(K.mat)[1])
    c.s = rep(dim(K.mat)[1],dim(K.mat)[2])
    while(max(abs(diag(u) %*% K.mat %*% v-r.s))>0.02 | max(abs( t(u)%*% K.mat %*% diag(v)-c.s))>0.02){
      v <- c(c.s/(t(u) %*% K.mat))
      u <- c(r.s/(K.mat %*% v))
    }
    
    K.mat.OT = diag(u) %*% K.mat %*% diag(v)
    
    K.mat.all = exp(-as.matrix(dist.mat)^2/quantile(dist.mat,0.5)^2)
    
    K.mat1=exp(-as.matrix(dist.mat)[1:dim(data.Y1)[2],1:dim(data.Y1)[2]]^2/quantile(dist.mat[1:dim(data.Y1)[2],1:dim(data.Y1)[2]],0.5)^2)
    K.mat2=exp(-as.matrix(dist.mat)[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]^2/quantile(dist.mat[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])],0.5)^2)
    
    
    
    #par(mfcol = c(1, 4), mar = c(4,2,1,1) + 0.1)
    
    
    #PCA
    sepl.svd = svds(data.Y2,k=r+1)
    out.sepl2 = sepl.svd$v[,1:r] 
    sepl.svd = svds(data.Y1,k=r+1)
    out.sepl1 = sepl.svd$v[,1:r] 
    
    #data.c=rbind(out.sepl1, out.sepl2)
    #plot(x=data.c[,1], y=data.c[,2], col = (c(rep("red",n),rep("blue",n))),pch=20, xlab="")
    #plot(x=data.c[,2], y=data.c[,3], col = (c(rep("red",n),rep("blue",n))),pch=20)
    
    #KPCA
    sep.svd = eigs(K.mat2,k=r+1)
    out.sep2 = sep.svd$vectors[,1:r+1] 
    sep.svd = eigs(K.mat1,k=r+1)
    out.sep1 = sep.svd$vectors[,1:r+1] 
    
    #data.c=rbind(out.sep1, out.sep2)
    #plot(x=data.c[,1], y=data.c[,2], col = (c(rep("red",n),rep("blue",n))),pch=20)
    #plot(x=data.c[,2], y=data.c[,3], col = (c(rep("red",n),rep("blue",n))),pch=20)
    
    #joint PCA
    pca.svd = svds(cbind(data.Y1,data.Y2), k=r)
    X.combined =pca.svd$v[,1:r] 
    out.pca2 = X.combined[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
    out.pca1 = X.combined[(1):(dim(data.Y1)[2]),]
    
    #data.c=rbind(out.pca1, out.pca2)
    #plot(x=data.c[,1], y=data.c[,2], col = (c(rep("red",n),rep("blue",n))),pch=20)
    #plot(x=data.c[,2], y=data.c[,3], col = (c(rep("red",n),rep("blue",n))),pch=20)
    
    #joint kPCA
    kpca.svd = eigs(K.mat.all, k=r+1)
    X.combined = kpca.svd$vectors[,1:r+1] 
    out.kpca2 = X.combined[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
    out.kpca1 = X.combined[(1):(dim(data.Y1)[2]),]
    
    #data.c=rbind(out.kpca1, out.kpca2)
    #plot(x=data.c[,1], y=data.c[,2], col = (c(rep("red",n),rep("blue",n))),pch=20)
    #plot(x=data.c[,2], y=data.c[,3], col = (c(rep("red",n),rep("blue",n))),pch=20)
    
    
    #CCA
    cca.svd = svds(t(data.Y1) %*% (data.Y2), k=r+1)
    out.cca2 = cca.svd$v[,1:r] 
    out.cca1 = cca.svd$u[,1:r]
    
    #data.c=rbind(out.cca1, out.cca2)
    #plot(x=data.c[,1], y=data.c[,2], col = (c(rep("red",n),rep("blue",n))),pch=20)
    #plot(x=data.c[,2], y=data.c[,3], col = (c(rep("red",n),rep("blue",n))),pch=20)  
    
    #LBDM
    alpha=1
    A = K.mat.rl
    D1 = rowSums(K.mat.rl)
    D2 = colSums(K.mat.rl)
    A = diag(D1^(-1/2)) %*% A %*% diag(D2^(-1/2))
    lbdm.svd = svds(A, k=r+1)
    X.combined = rbind(lbdm.svd$u[,1:r+1], lbdm.svd$v[,1:r+1])
    V = diag(c(D1,D2)^(-1/2)) %*% X.combined %*% diag(lbdm.svd$d[1:r+1]^alpha)
    out.lbdm2 = V[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
    out.lbdm1 = V[(1):(dim(data.Y1)[2]),]
    
    #data.c=rbind(out.lbdm1, out.lbdm2)
    #plot(x=data.c[,1], y=data.c[,2], col = (c(rep("red",n),rep("blue",n))),pch=20)
    #plot(x=data.c[,2], y=data.c[,3], col = (c(rep("red",n),rep("blue",n))),pch=20)
    
    #roseland
    
    D.r = rowSums(t(K.mat.rl) %*% (K.mat.rl))
    r.svd.r = svds(diag((D.r)^(-1/2)) %*% t(K.mat.rl), k=r+1)
    out.rl2=diag((D.r)^(-1/2)) %*% r.svd.r$u[,2:(r+1)] %*% diag(r.svd.r$d[2:(r+1)])^2
    D.r = rowSums((K.mat.rl) %*% t(K.mat.rl))
    r.svd.r = svds(diag((D.r)^(-1/2)) %*% (K.mat.rl), k=r+1)
    out.rl1=diag((D.r)^(-1/2)) %*% r.svd.r$u[,2:(r+1)] %*% diag(r.svd.r$d[2:(r+1)])^2
    
    #data.c=rbind(out.rl1, out.rl2)
    #plot(x=data.c[,1], y=data.c[,2], col = (c(rep("red",n),rep("blue",n))),pch=20)
    #plot(x=data.c[,2], y=data.c[,3], col = (c(rep("red",n),rep("blue",n))),pch=20)
    
    #Prop-0
    prop.svd = svds(K.mat.OT, k=r+2)
    out.prop2 = sqrt(dim(data.X2)[2])*prop.svd$v[,1:r+1]
    out.prop1 = sqrt(dim(data.X1)[2])*prop.svd$u[,1:r+1]
    
    #data.c=rbind(out.prop1, out.prop2)
    #plot(x=data.c[,1], y=data.c[,2], col = (c(rep("red",n),rep("blue",n))),pch=20)
    #plot(x=data.c[,2], y=data.c[,3], col = (c(rep("red",n),rep("blue",n))),pch=20)
    
    #Prop-1
    prop.svd = svds(K.mat.OT, k=r+2)
    out.prop22 = sqrt(dim(data.X2)[2])*prop.svd$v[,1:r+1] %*% diag(prop.svd$d[2:(r+1)]^1)
    out.prop12 = sqrt(dim(data.X1)[2])*prop.svd$u[,1:r+1] %*% diag(prop.svd$d[2:(r+1)]^1)
    
    #data.c=rbind(out.prop12, out.prop22)
    #plot(x=data.c[,1], y=data.c[,2], col = (c(rep("red",n),rep("blue",n))),pch=20)
    #plot(x=data.c[,2], y=data.c[,3], col = (c(rep("red",n),rep("blue",n))),pch=20)
    
    
    error12[NN,]=c(index.DB(rbind(out.sepl1,out.sepl2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.sep1,out.sep2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.pca1,out.pca2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.kpca1,out.kpca2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.cca1,out.cca2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.lbdm1,out.lbdm2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.rl1,out.rl2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.prop1,out.prop2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.prop12,out.prop22), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB)
    
    error1[NN,]=c(knn_dist(rbind(data.X1,data.X2),rbind(out.sepl1,out.sepl2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.sep1,out.sep2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.pca1,out.pca2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.kpca1,out.kpca2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.cca1,out.cca2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.lbdm1,out.lbdm2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.rl1,out.rl2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.prop1,out.prop2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.prop12,out.prop22)))
    
    print(NN)
  }
  print(pp)
  error12.out[pp,] = colMeans(error12,na.rm = T)
  error1.out[pp,] = colMeans(error1,na.rm = T)
}


data_p = data.frame(j.index = c(error1.out),
                    method = factor(rep(c("pca", "kpca", "j-pca", "j-kpca","seurat",
                                          "lbdm","rl","EOT-0","EOT-1"), 
                                        each=length(sig.v))),
                    sigma = rep(signif(sig.v,3), times = 9))

data_p$method <- factor(data_p$method, levels = c("EOT-0","EOT-1","lbdm","rl","seurat","j-pca", "j-kpca", "kpca", "pca"))

plt <- ggplot(data=data_p, aes(x=sigma, y=j.index), group=factor(method)) +
  geom_line(aes(linetype=method,color=method))+
  geom_point(aes(shape=method,color=method) ,size=5) + ylab("concordance (Jaccard index)") + xlab("tau")+
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #600*500

plt 

data_p = data.frame(rand.index = log(c(error12.out)),
                    method = factor(rep(c("pca", "kpca", "j-pca", "j-kpca","seurat",
                                          "lbdm","rl","EOT-0","EOT-3"), 
                                        each=length(sig.v))),
                    sigma = rep(sig.v, times = 9))

data_p$method <- factor(data_p$method, levels = c("EOT-0","EOT-1","lbdm","rl","seurat","j-pca", "j-kpca", "kpca", "pca"))

plt <- ggplot(data=data_p, aes(x=sigma, y=rand.index), group=factor(method)) +
  geom_line(aes(linetype=method,color=method))+
  geom_point(aes(shape=method,color=method) ,size=5) + ylab("alignment (D-B index)") + xlab("tau")+
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #600*500

plt 




### Setting 2


sig.v=(seq(0.001,1, length.out= 20)) 
time=matrix(ncol=9,nrow=500)
error12=matrix(ncol=9,nrow=500)
error12.out=matrix(ncol=9,nrow=length(sig.v))
error2=matrix(ncol=9,nrow=500)
error2.out=matrix(ncol=9,nrow=length(sig.v))
error1=matrix(ncol=9,nrow=500)
error1.out=matrix(ncol=9,nrow=length(sig.v))
n=600
theta=n^(2/5)

for(pp in 1:length(sig.v)){
  p=1000
  sig1 = theta*0.05
  sigma=sig.v[pp]
  for(NN in 1:5){
    
    data.X1=data[sample(dim(data)[1],n),]
    data.X2=data[sample(dim(data)[1],n),]
    data.Z=matrix(rnorm(n*p,0,sig1), ncol=p)
    data.Y1 = data.Z
    data.Y1[,1:r] = data.Y1[,1:r]+as.matrix(data.X1*theta)
    data.Y1[,1] = data.Y1[,1] + theta*3
    data.Y1 = data.Y1*3
    
    data.Z=matrix(rnorm(n*p,0,sig1), ncol=p)
    data.Y2 = data.Z
    data.Y2[,1:r] = data.Y2[,1:r]+as.matrix(data.X2*theta)
    data.Y2[1:floor(n/3),2:r] = data.Y2[1:floor(n/3),2:r]+matrix(rnorm(floor(n/3)*2,0,sig1*3), ncol=2)
    data.Y2[(floor(n/3)+1):floor(2*n/3),1:r] = data.Y2[(floor(n/3)+1):floor(2*n/3),1:r]+matrix(rnorm(floor(n/3)*r,0,2*sig1), ncol=r)
    data.Y2[,(r+1):(p)] = data.Y2[,(r+1):(p)]+matrix(runif((n)*length((r+1):(p)),sigma*theta/2,sigma*theta),ncol=length((r+1):(p)))
    
    #various kernel matrices
    data.Y1=t(data.Y1)
    data.Y2=t(data.Y2)
    dist.mat = Dist(t(cbind(data.Y1, data.Y2)))
    cross.dist = as.matrix(dist.mat)[1:dim(data.Y1)[2], (dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]
    
    K.mat = exp(-cross.dist^2/quantile(cross.dist,0.5)^2) 
    K.mat.rl =  K.mat
    
    ##Sinkhorn-Knopp algorithm
    u=rep(1,dim(K.mat)[1])
    v=rep(1,dim(K.mat)[2])
    r.s = rep(dim(K.mat)[2],dim(K.mat)[1])
    c.s = rep(dim(K.mat)[1],dim(K.mat)[2])
    while(max(abs(diag(u) %*% K.mat %*% v-r.s))>0.02 | max(abs( t(u)%*% K.mat %*% diag(v)-c.s))>0.02){
      v <- c(c.s/(t(u) %*% K.mat))
      u <- c(r.s/(K.mat %*% v))
    }
    
    K.mat.OT = diag(u) %*% K.mat %*% diag(v)
    
    K.mat.all = exp(-as.matrix(dist.mat)^2/quantile(dist.mat,0.5)^2)
    
    K.mat1=exp(-as.matrix(dist.mat)[1:dim(data.Y1)[2],1:dim(data.Y1)[2]]^2/quantile(dist.mat[1:dim(data.Y1)[2],1:dim(data.Y1)[2]],0.5)^2)
    K.mat2=exp(-as.matrix(dist.mat)[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])]^2/quantile(dist.mat[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2])],0.5)^2)
    
    
    
    #par(mfcol = c(1, 3))
    
    
    #PCA
    sepl.svd = svds(data.Y2,k=r+1)
    out.sepl2 = sepl.svd$v[,1:r] 
    sepl.svd = svds(data.Y1,k=r+1)
    out.sepl1 = sepl.svd$v[,1:r] 
    
    #data.c=rbind(out.sepl1, out.sepl2)
    #scatterplot3d(x=data.c[,1], y=data.c[,2], z=data.c[,3], angle=45, color = (c(rep("red",n),rep("blue",n))),pch=20,
    #              xlab=NA, ylab=NA, zlab=NA, label.tick.marks = FALSE)
    
    #KPCA
    sep.svd = eigs(K.mat2,k=r+1)
    out.sep2 = sep.svd$vectors[,1:r+1] 
    sep.svd = eigs(K.mat1,k=r+1)
    out.sep1 = sep.svd$vectors[,1:r+1] 
    
    #data.c=rbind(out.sep1, out.sep2)
    #scatterplot3d(x=data.c[,1], y=data.c[,2], z=data.c[,3], angle=45, color = (c(rep("red",n),rep("blue",n))),pch=20,
    #              xlab=NA, ylab=NA, zlab=NA, label.tick.marks = FALSE)
    
    #joint PCA
    pca.svd = svds(cbind(data.Y1,data.Y2), k=r)
    X.combined =pca.svd$v[,1:r] 
    out.pca2 = X.combined[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
    out.pca1 = X.combined[(1):(dim(data.Y1)[2]),]
    
    #data.c=rbind(out.pca1, out.pca2)
    #scatterplot3d(x=data.c[,1], y=data.c[,2], z=data.c[,3], angle=45, color = (c(rep("red",n),rep("blue",n))),pch=20,
    #              xlab=NA, ylab=NA, zlab=NA, label.tick.marks = FALSE)
    
    #joint kPCA
    kpca.svd = eigs(K.mat.all, k=r+1)
    X.combined = kpca.svd$vectors[,1:r+1] 
    out.kpca2 = X.combined[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
    out.kpca1 = X.combined[(1):(dim(data.Y1)[2]),]
    
    #data.c=rbind(out.kpca1, out.kpca2)
    #scatterplot3d(x=data.c[,1], y=data.c[,2], z=data.c[,3], angle=45, color = (c(rep("red",n),rep("blue",n))),pch=20,
    #              xlab=NA, ylab=NA, zlab=NA, label.tick.marks = FALSE)
    
    
    #CCA
    cca.svd = svds(t(data.Y1) %*% (data.Y2), k=r+1)
    out.cca2 = cca.svd$v[,1:r] 
    out.cca1 = cca.svd$u[,1:r]
    
    #data.c=rbind(out.cca1, out.cca2)
    #scatterplot3d(x=data.c[,1], y=data.c[,2], z=data.c[,3], angle=45, color = (c(rep("red",n),rep("blue",n))),pch=20,
    #              xlab=NA, ylab=NA, zlab=NA, label.tick.marks = FALSE)
    
    #LBDM
    alpha=1
    A = K.mat.rl
    D1 = rowSums(K.mat.rl)
    D2 = colSums(K.mat.rl)
    A = diag(D1^(-1/2)) %*% A %*% diag(D2^(-1/2))
    lbdm.svd = svds(A, k=r+1)
    X.combined = rbind(lbdm.svd$u[,1:r+1], lbdm.svd$v[,1:r+1])
    V = diag(c(D1,D2)^(-1/2)) %*% X.combined %*% diag(lbdm.svd$d[1:r+1]^alpha)
    out.lbdm2 = V[(dim(data.Y1)[2]+1):(dim(data.Y1)[2]+dim(data.Y2)[2]),]
    out.lbdm1 = V[(1):(dim(data.Y1)[2]),]
    
    #data.c=rbind(out.lbdm1, out.lbdm2)
    #scatterplot3d(x=data.c[,1], y=data.c[,2], z=data.c[,3], angle=45, color = (c(rep("red",n),rep("blue",n))),pch=20,
    #              xlab=NA, ylab=NA, zlab=NA, label.tick.marks = FALSE)
    
    #roseland
    
    D.r = rowSums(t(K.mat.rl) %*% (K.mat.rl))
    r.svd.r = svds(diag((D.r)^(-1/2)) %*% t(K.mat.rl), k=r+1)
    out.rl2=diag((D.r)^(-1/2)) %*% r.svd.r$u[,2:(r+1)] %*% diag(r.svd.r$d[2:(r+1)])^2
    D.r = rowSums((K.mat.rl) %*% t(K.mat.rl))
    r.svd.r = svds(diag((D.r)^(-1/2)) %*% (K.mat.rl), k=r+1)
    out.rl1=diag((D.r)^(-1/2)) %*% r.svd.r$u[,2:(r+1)] %*% diag(r.svd.r$d[2:(r+1)])^2
    
    #data.c=rbind(out.rl1, out.rl2)
    #scatterplot3d(x=data.c[,1], y=data.c[,2], z=data.c[,3], angle=45, color = (c(rep("red",n),rep("blue",n))),pch=20,
    #              xlab=NA, ylab=NA, zlab=NA, label.tick.marks = FALSE)
    
    #Prop-0
    prop.svd = svds(K.mat.OT, k=r+2)
    out.prop2 = sqrt(dim(data.X2)[2])*prop.svd$v[,1:r+1]
    out.prop1 = sqrt(dim(data.X1)[2])*prop.svd$u[,1:r+1]
    
    #data.c=rbind(out.prop1, out.prop2)
    #scatterplot3d(x=data.c[,1], y=data.c[,2], z=data.c[,3], angle=45, color = (c(rep("red",n),rep("blue",n))),pch=20,
    #              xlab=NA, ylab=NA, zlab=NA, label.tick.marks = FALSE)
    
    #Prop-1
    prop.svd = svds(K.mat.OT, k=r+2)
    out.prop22 = sqrt(dim(data.X2)[2])*prop.svd$v[,1:r+1] %*% diag(prop.svd$d[2:(r+1)]^1)
    out.prop12 = sqrt(dim(data.X1)[2])*prop.svd$u[,1:r+1] %*% diag(prop.svd$d[2:(r+1)]^1)
    
    #data.c=rbind(out.prop12, out.prop22)
    #scatterplot3d(x=data.c[,1], y=data.c[,2], z=data.c[,3], angle=45, color = (c(rep("red",n),rep("blue",n))),pch=20,
    #              xlab=NA, ylab=NA, zlab=NA, label.tick.marks = FALSE)
    
    
    error12[NN,]=c(index.DB(rbind(out.sepl1,out.sepl2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.sep1,out.sep2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.pca1,out.pca2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.kpca1,out.kpca2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.cca1,out.cca2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.lbdm1,out.lbdm2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.rl1,out.rl2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.prop1,out.prop2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.prop12,out.prop22), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB)
    
    error1[NN,]=c(knn_dist(rbind(data.X1,data.X2),rbind(out.sepl1,out.sepl2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.sep1,out.sep2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.pca1,out.pca2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.kpca1,out.kpca2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.cca1,out.cca2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.lbdm1,out.lbdm2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.rl1,out.rl2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.prop1,out.prop2)),
                  knn_dist(rbind(data.X1,data.X2),rbind(out.prop12,out.prop22)))
    
    print(NN)
  }
  print(pp)
  error12.out[pp,] = colMeans(error12,na.rm = T)
  error1.out[pp,] = colMeans(error1,na.rm = T)
}


data_p = data.frame(j.index = c(error1.out),
                    method = factor(rep(c("pca", "kpca", "j-pca", "j-kpca","seurat",
                                          "lbdm","rl","EOT-0","EOT-1"), 
                                        each=length(sig.v))),
                    sigma = rep(signif(sig.v,3), times = 9))

data_p$method <- factor(data_p$method, levels = c("EOT-0","EOT-1","lbdm","rl","seurat","j-pca", "j-kpca", "kpca", "pca"))

plt <- ggplot(data=data_p, aes(x=sigma, y=j.index), group=factor(method)) +
  geom_line(aes(linetype=method,color=method))+
  geom_point(aes(shape=method,color=method) ,size=5) + ylab("concordance (Jaccard index)") + xlab("gamma")+
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #600*500

plt 

data_p = data.frame(rand.index = log(c(error12.out)),
                    method = factor(rep(c("pca", "kpca", "j-pca", "j-kpca","seurat",
                                          "lbdm","rl","EOT-0","EOT-1"), 
                                        each=length(sig.v))),
                    sigma = rep(sig.v, times = 9))

data_p$method <- factor(data_p$method, levels = c("EOT-0","EOT-1","lbdm","rl","seurat","j-pca", "j-kpca", "kpca", "pca"))

plt <- ggplot(data=data_p, aes(x=sigma, y=rand.index), group=factor(method)) +
  geom_line(aes(linetype=method,color=method))+
  geom_point(aes(shape=method,color=method) ,size=5) + ylab("alignment (D-B index)") + xlab("gamma")+
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #600*500

plt 
