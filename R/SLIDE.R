# This function is to estimate pairwise linkage disequilibrium with EM algorithm

require(MASS)

SILDE<-function(pheno_geno,method,LD=NULL,num_per=200)
{

  n_gv<-ncol(pheno_geno)-1

  fa_case<-pheno_geno[which(pheno_geno[,1]==1),-1]

  n_case<-nrow(fa_case)

  fa_con<-pheno_geno[which(pheno_geno[,1]==0),-1]
  n_con<-nrow(fa_con)

  U<-(colSums(fa_case)/n_case-colSums(fa_con)/n_con)/(1/n_case+1/n_con)


  if(method=="Cov") V1<-cov(fa_con)   else V1<-2*ld.em(fa_con)[[1]]

  LD_matrix<-list()

 LD_matrix[[1]]<-V1*n_case*n_con/(n_case+n_con)

  if(length(LD)>0) LD_matrix[[2]]<-2*LD*n_case*n_con/(n_case+n_con)


 LD_inverse<-list()

 TT<-NULL

 for(ii in 1:length(LD_matrix))
 {
   LD_inverse[[ii]]<-ginv(LD_matrix[[ii]])

   TT<-c(TT,U%*%LD_inverse[[ii]]%*%U)

 }

  for(B in 1:num_per)
  {

    fa_per<-permuC(pheno_geno[,1],pheno_geno[,-1])

    fa_case<-fa_per[which(fa_per[,1]==1),-1]

    fa_con<-fa_per[which(fa_per[,1]==0),-1]

    U<-(colSums(fa_case)/n_case-colSums(fa_con)/n_con)/(1/n_case+1/n_con)

    tt0<-NULL

    for(ii in 1:length(LD_inverse)) tt0<-c(tt0,U%*%LD_inverse[[ii]]%*%U)


    TT<-rbind(TT,tt0)
  }

  all<-(num_per+1-apply(TT,2,rank))/num_per

  if(length(LD)>0) p_asy_out_cov<-1-pchisq(TT[1,2],n_gv) else p_asy_out_cov<-NULL

   p_asy_cov<-1-pchisq(TT[1,1],n_gv)

  p.value<-c(all[1,],p_asy_cov,p_asy_out_cov)

  if(length(p.value)>2) {

    SLIDE.tem<-c(TT[1,1],p.value[c(1,3)]); names(SLIDE.tem)<-c("SLIDE","p.value.per","p.value.asy")
    SLIDE.pop<-c(TT[1,2],p.value[c(1,4)]); names(SLIDE.pop)<-c("SLIDE.pop","p.value.pop.per","p.value.pop.asy")

    statis<-list(SLIDE.tem,SLIDE.pop)
  } else {
    SLIDE.tem<-c(TT[1,1],p.value); names(SLIDE.tem)<-c("SLIDE","p.value.per","p.value.asy")
    statis<-list(SLIDE.tem)
  }


  return(statis)
}
