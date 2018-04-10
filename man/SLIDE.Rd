\name{SLIDE}
\title{Retrospective Score Test for Case-Control Association Study}
\description{
To explain more heritability in GWAS, the score test statistics incorporating linkage disequilibrium information in the retrospective perspective for case-control studies are proposed. To be specific, the score is defined as the difference of the average multi-locus genotypes between cases and controls and then its  variance-covariance matrix involves linkage disequilibrium information. The essential difference of the variance-covariance matrix of the proposed test with Score test exists, despite that the forms are similar. A noticeable merit/feature of SLDE is that it could borrow the strength from a database with several thousands to hundreds of thousands size to improve the power for detecting association.}

\usage{SILDE(pheno_geno,method,LD=NULL,num_per=200)}
\arguments{
\item{pheno_geno}{  Matrix in which the first column is the phenotype (0=control, 1=case) and others are the corresponding genotype (0, 1, or 2). Each row represents an individual.}

\item{LD}{ The known LD of the underling population of the cases and controls. By default, argument LD=NULL meaning that the LD of the underling population is unknown.}

\item{method}{ The method is to calculate the LD using sample data. Argument method="Cov" meaning that the variance-covariance matrix is calculated by function cov(); meanwhile, method="EM" meaning that the variance-covariance matrix is calculated by EM algorithm.}

\item{num_per}{ Positive integer indicating the number of permutations (200 by default).}
}
\details{
The results with argument method="Cov" are similar to the results with method="EM", whereas the calculation with method="Cov" is simple and fast.

The asymptotical p-value of SILDE may be little inflated when the sample size is not enough big, especially for rare variants.}
\value{ A vector with the following elements:

\item{SILDE}{ Statistic SILDE with the estimated LD using controls, the corresponding asymptotical p-value and permuted p-value.}

\item{SILDE.pop}{ Statistic SILDE with the known LD of underling population, the corresponding asymptotical p-value and permuted p-value. If LD=NULL, then this is NULL.}

}

\references{Chan Wang, Leiming Sun, Liming Li and Yue-Qing Hu, Retrospective Score Test for Case-Control Association Study}
\author{ Chan Wang and Yue-Qing Hu }
\examples{

library(MASS)

set.seed(1234)

genotype<-matrix(sample(c(0,1,2),5000,replace=TRUE),500,10) ### 500 individuals and 10 SNPs

phenotype<-c(rep(1,200),rep(0,300)) ### 200 cases and 300 controls

Statistics<-SILDE(cbind(phenotype,genotype),method=c("Cov"),LD=NULL,num_per=200)

SLIDE<-Statistics[[1]]

## 21.2053, 0.0200, 0.0197

}

