# This function is to estimate pairwise linkage disequilibrium with EM algorithm

ld.em<-function(genotype)
{

  n_site<-ncol(genotype)

  dat<-genotype
  n<-nrow(dat)

  LD0<-LD1<-LD2<-diag((colSums(dat)/2/n)*(1-(colSums(dat)/2/n)))
  for(i in 1:(n_site-1))
  {
    for(j in c(i:n_site)[which(c(i:n_site)!=i)])
    {
      p_A<-(2*n-sum(dat[,i]))/2/n
      p_B<-(2*n-sum(dat[,j]))/2/n
      p_AB<-p_A*p_B
      p_Ab<-p_A*(1-p_B)
      p_aB<-(1-p_A)*p_B
      p_ab<-(1-p_A)*(1-p_B)

      n_AABB<-length(which(dat[,i]==0 & dat[,j]==0))
      n_AABb<-length(which(dat[,i]==0 & dat[,j]==1))
      n_AaBB<-length(which(dat[,i]==1 & dat[,j]==0))
      n_AaBb<-length(which(dat[,i]==1 & dat[,j]==1))
      n_AAbb<-length(which(dat[,i]==0 & dat[,j]==2))
      n_Aabb<-length(which(dat[,i]==1 & dat[,j]==2))
      n_aabb<-length(which(dat[,i]==2 & dat[,j]==2))
      n_aaBb<-length(which(dat[,i]==2 & dat[,j]==1))
      n_aaBB<-length(which(dat[,i]==2 & dat[,j]==0))

      p_AB_i<-p_AB
      p_ab_i<-p_ab
      p_Ab_i<-p_Ab
      p_aB_i<-p_aB

      for(t in 1:100)
      {
        a<-p_AB_i[t]*p_ab_i[t]/(p_AB_i[t]*p_ab_i[t]+p_Ab_i[t]*p_aB_i[t])

        b<-p_Ab_i[t]*p_aB_i[t]/(p_AB_i[t]*p_ab_i[t]+p_Ab_i[t]*p_aB_i[t])

        p_AB_i<-rbind(p_AB_i,(2*n_AABB+n_AABb+n_AaBB+n_AaBb*a)/2/n)

        p_ab_i<-rbind(p_ab_i,(2*n_aabb+n_aaBb+n_Aabb+n_AaBb*a)/2/n)

        p_Ab_i<-rbind(p_Ab_i,(2*n_AAbb+n_AABb+n_Aabb+n_AaBb*b)/2/n)

        p_aB_i<-rbind(p_aB_i,(2*n_aaBB+n_aaBb+n_AaBB+n_AaBb*b)/2/n)

        aa<-max(abs(c(p_AB_i[t+1]-p_AB_i[t],p_Ab_i[t+1]-p_Ab_i[t],p_aB_i[t+1]-p_aB_i[t],p_ab_i[t+1]-p_ab_i[t])))

        if(aa<=0.00001)   break
      }
      LD0[i,j]<-p_AB_i[t+1]*p_ab_i[t+1]-p_Ab_i[t+1]*p_aB_i[t+1]
      LD0[j,i]<-LD0[i,j]

      LD1[i,j]<-LD1[j,i]<-abs(LD0[i,j])/(min(p_A*(1-p_B),(1-p_A)*p_B)*sign(LD0[i,j]>0)+min((1-p_A)*(1-p_B), p_A*p_B)*sign(LD0[i,j]<0)+sign(LD0[i,j]==0))

      LD2[i,j]<- LD2[j,i]<-LD0[i,j]^2/(p_A*(1-p_B)*(1-p_A)*p_B)

    }
  }
list(LD0,LD1,LD2)
}


