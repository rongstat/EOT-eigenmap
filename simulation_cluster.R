######## simulation for KSBE: denoising (numerical results)

library(Rfast)
library(RSpectra)
library(BiocNeighbors)
library(fossil)
library(ggplot2)
library(cluster)
library(uwot)
###
r=6
c.prop = c(1,1,1,1,1,1)
c.prop = c.prop/sum(c.prop)#class proportion
rho=5

#set mean vectors
mu=rep(0,r)
mu1=mu
mu1[1]=rho
mu2=mu
mu2[2]=rho
mu3=mu
mu3[3]=rho
mu4=mu
mu4[4]=rho
mu5=mu
mu5[5]=rho
mu6=mu
mu6[6]=rho
sigma=diag(rep(1,r))


error=matrix(ncol=9,nrow=500)
error12=matrix(ncol=9,nrow=500)
sig.v=seq(1,3, length.out= 20)
error.out=matrix(ncol=9,nrow=length(sig.v))
error12.out=matrix(ncol=9,nrow=length(sig.v))

for(pp in 1:length(sig.v)){
  theta=sig.v[pp]
  n0=600
  n=sum(round(n0*c.prop))
  p=1000
  #epsilon = (log(n)/n)^(1/(4*r+5))
  for(NN in 1:20){
    
    
    #setting 1:
    #generate 6 clusters
    data1=rmvnorm(round(n0*c.prop[1]),mu1,sigma)
    data2=rmvnorm(round(n0*c.prop[2]),mu2,sigma)
    data3=rmvnorm(round(n0*c.prop[3]),mu3,sigma)
    data4=rmvnorm(round(n0*c.prop[4]),mu4,sigma)
    data5=rmvnorm(round(n0*c.prop[5]),mu5,sigma)
    data6=rmvnorm(round(n0*c.prop[6]),mu6,sigma)
    # data set X
    data.X1=rbind(data1,data2,data3,data4,data5,data6)
    data.Z1=rmvnorm(n,rep(0,p),diag(rep(0.5,p)))
    data.Y1 = data.Z1
    data.Y1[,1:6] = data.Y1[,1:6]+data.X1*theta
    data.Y1[,3] = data.Y1[,3] + theta*5
    data.Y1[,2] = data.Y1[,2] + theta*5
    #data.Y1=5*data.Y1
    info1=rep(1:6, times = round(n0*c.prop))
    #plot(umap(data.Y1,spread=2),col=info1)
    
    data1=rmvnorm(round(n0*c.prop[1]),mu1,sigma)
    data2=rmvnorm(round(n0*c.prop[2]),mu2,sigma)
    data3=rmvnorm(round(n0*c.prop[3]),mu3,sigma)
    data4=rmvnorm(round(n0*c.prop[4]),mu4,sigma)
    data5=rmvnorm(round(n0*c.prop[5]),mu5,sigma)
    data6=rmvnorm(round(n0*c.prop[6]),mu6,sigma)
    # data set X
    data.X2=rbind(data1,data2,data3,data4,data5,data6)
    data.Z2=rmvnorm(n,rep(0,p),diag(rep(1,p)))
    data.Y2 = data.Z2
    data.Y2[,1:6] = data.Y2[,1:6]+data.X2*theta
    data.Y2[1:floor(n/3),2:r] = data.Y2[1:floor(n/3),2:r]+matrix(rnorm(floor(n/3)*5,0,3), ncol=5)
    data.Y2[(floor(n/3)+1):floor(2*n/3),1:r] = data.Y2[(floor(n/3)+1):floor(2*n/3),1:r]+matrix(rnorm(floor(n/3)*r,0,2), ncol=r)
    data.Y2[,(r+1):(p)] = data.Y2[,(r+1):(p)]+matrix(runif((n)*length((r+1):(p)),theta/2,theta),ncol=length((r+1):(p)))
    
    #data.Y2[,(r+1):p] = data.Y2[,(r+1):p]+matrix(runif(n*(p-r),-10,10),ncol=p-r) #noise in Y2
    info2=rep(1:6, times =  round(n0*c.prop))
    
    
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
    
    
    
    #par(mfrow = c(2, 4), mar = c(4,1,1,1) + 0.1)
    
    #PCA
    sepl.svd = svds(data.Y2,k=r+1)
    out.sepl2 = sepl.svd$v[,1:r] %*% diag(sepl.svd$d[1:r])^(1/2)
    sepl.svd = svds(data.Y1,k=r+1)
    out.sepl1 = sepl.svd$v[,1:r] %*% diag(sepl.svd$d[1:r])^(1/2)
    cls=kmeans(rbind(out.sepl1,out.sepl2),centers=6)$cluster
    error[NN,1] = rand.index(cls,c(info1,info2))
    
    
    #KPCA
    sep.svd = eigs(K.mat2,k=r+1)
    out.sep2 = sep.svd$vectors[,1:r+1] %*% diag(sep.svd$values[1:r+1])^(1/2)
    sep.svd = eigs(K.mat1,k=r+1)
    out.sep1 = sep.svd$vectors[,1:r+1] %*% diag(sep.svd$values[1:r+1])^(1/2)
    cls=kmeans(rbind(out.sep1,out.sep2),centers=6)$cluster
    error[NN,2] = rand.index(cls,c(info1,info2))
    
    
    #joint PCA
    pca.svd = svds(cbind(data.Y1,data.Y2), k=r)
    X.jpca =pca.svd$v[,1:r] %*% diag(pca.svd$d[1:r])^(1/2)
    cls=kmeans(X.jpca,centers=6)$cluster
    error[NN,3] = rand.index(cls,c(info1,info2))
    
    #joint kPCA
    kpca.svd = eigs(K.mat.all, k=r+1)
    X.jkpca = kpca.svd$vectors[,1:r+1] %*% diag(kpca.svd$values[1:r+1])^(1/2)
    cls=kmeans(X.jkpca,centers=6)$cluster
    error[NN,4] = rand.index(cls,c(info1,info2))
    
    #Seurat
    cca.svd = svds(t(data.Y1) %*% (data.Y2), k=r+1)
    out.cca2 = cca.svd$v[,1:r] 
    out.cca1 = cca.svd$u[,1:r]
    X.combined = rbind(out.cca1,out.cca2)
    cls=kmeans(X.combined,centers=6)$cluster
    error[NN,5] = rand.index(cls,c(info1,info2))
    
    
    #LBDM
    alpha=1
    A = K.mat.rl
    D1 = rowSums(K.mat.rl)
    D2 = colSums(K.mat.rl)
    A = diag(D1^(-1/2)) %*% A %*% diag(D2^(-1/2))
    lbdm.svd = svds(A, k=r+1)
    X.combined = rbind(lbdm.svd$u[,1:r+1], lbdm.svd$v[,1:r+1])
    V = diag(c(D1,D2)^(-1/2)) %*% X.combined %*% diag(lbdm.svd$d[1:r+1]^alpha)
    out.lbdm = V
    cls=kmeans(V,centers=6)$cluster
    error[NN,6] = rand.index(cls,c(info1,info2))
    
    #roseland
    D.r = rowSums(t(K.mat.rl) %*% (K.mat.rl))
    r.svd.r = svds(diag((D.r)^(-1/2)) %*% t(K.mat.rl), k=r+1)
    out.rl2=diag((D.r)^(-1/2)) %*% r.svd.r$u[,2:(r+1)] %*% diag(r.svd.r$d[2:(r+1)])^2
    D.r = rowSums((K.mat.rl) %*% t(K.mat.rl))
    r.svd.r = svds(diag((D.r)^(-1/2)) %*% (K.mat.rl), k=r+1)
    out.rl1=diag((D.r)^(-1/2)) %*% r.svd.r$u[,2:(r+1)] %*% diag(r.svd.r$d[2:(r+1)])^2
    cls=kmeans(rbind(out.rl1,out.rl2),centers=6)$cluster
    error[NN,7] = rand.index(cls,c(info1,info2))
    
    #Prop
    #Prop-0
    prop.svd = svds(K.mat.OT, k=r+2)
    out.prop2 = sqrt(dim(data.X2)[2])*prop.svd$v[,1:r+1]
    out.prop1 = sqrt(dim(data.X1)[2])*prop.svd$u[,1:r+1]
    cls=kmeans(rbind(out.prop1,out.prop2),centers=6)$cluster
    error[NN,8] = rand.index(cls,c(info1,info2))
    
    #Prop-1
    prop.svd = svds(K.mat.OT, k=r+2)
    out.prop22 = sqrt(dim(data.X2)[2])*prop.svd$v[,1:r+1] %*% diag(prop.svd$d[2:(r+1)]^1)
    out.prop12 = sqrt(dim(data.X1)[2])*prop.svd$u[,1:r+1] %*% diag(prop.svd$d[2:(r+1)]^1)
    cls=kmeans(rbind(out.prop12,out.prop22),centers=6)$cluster
    error[NN,9] = rand.index(cls,c(info1,info2))
    
    error12[NN,]=c(index.DB(rbind(out.sepl1,out.sepl2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.sep1,out.sep2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(X.jpca, cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(X.jkpca, cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.cca1,out.cca2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB((out.lbdm), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.rl1,out.rl2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.prop1,out.prop2), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB,
                   index.DB(rbind(out.prop12,out.prop22), cl=as.numeric(c(rep(1,n), rep(2,n))))$DB)
    
    print(NN)
  }
  error.out[pp,] = colMeans(error,na.rm=T)
  error12.out[pp,] = colMeans(error12,na.rm = T)
}

#error.out

data = data.frame(rand.index = c(error.out),
                  method = factor(rep(c("pca", "kpca", "j-pca", "j-kpca","seurat",
                                        "lbdm","rl","EOT-0","EOT-1"), 
                                      each=length(sig.v))),
                  theta = rep(sig.v, times = 9))

data$method <- factor(data$method, levels = c("EOT-0","EOT-1","lbdm","rl","seurat","j-pca", "j-kpca", "kpca", "pca"))

ggplot(data, aes(x=theta, y=rand.index, group=factor(method), colour = method)) +
  geom_line(aes(linetype=method))+xlab("theta")+
  geom_point(aes(shape=method) ,size=5) + ylab("rand index") + 
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #600*500


data = data.frame(rand.index = c(error12.out),
                  method = factor(rep(c("pca", "kpca", "j-pca", "j-kpca","seurat",
                                        "lbdm","rl","EOT-0","EOT-1"), 
                                      each=length(sig.v))),
                  theta = rep(sig.v, times = 9))

data$method <- factor(data$method, levels = c("EOT-0","EOT-1","lbdm","rl","seurat","j-pca", "j-kpca", "kpca", "pca"))

ggplot(data, aes(x=theta, y=log(rand.index), group=factor(method), colour = method)) +
  geom_line(aes(linetype=method))+xlab("theta")+
  geom_point(aes(shape=method) ,size=5) + ylab("log(D-B index)") + 
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) #600*500


