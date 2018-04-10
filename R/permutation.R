## The permutation for case-control studies

permuC<-function(phenotype,genotype)
{
  n<-nrow(genotype)
  size1<-sum(phenotype==1)
  size2<-sum(phenotype==0)
  fa_per<-genotype[sample(c(1:n),n),]
  fa_per<-cbind(matrix(c(rep(1,size1),rep(0,size2)),ncol=1),fa_per)

  fa_per
}

