\name{LD-EM}
\title{Estimating LD using EM algorithm}
\description{
This function is to estimate pairwise linkage disequilibrium with EM algorithm. }

\usage{ld.em(genotype)}
\arguments{
\item{genotype}{ numeric matrix or data frame with genotype data coded as 0, 1, 2. No missing data is allowed.}
}
\details{There is no missing data.}
\value{The list contains the LD matrix, \eqn{D'} matrix and \eqn{r^2} matrix.}

\references{Chan Wang, Leiming Sun, Liming Li and Yue-Qing Hu, Linkage Disequilibrium on Improving Power in Case-Control Association Study}
\author{ Chan Wang and Yue-Qing Hu }
\examples{
genotype<-matrix(sample(c(0,1,2),5000,replace=TRUE),500,10) ### 500 individuals and 10 SNPs

LD<-ld.em(genotype)

LD_orignal<-LD[[1]]

LD_D<-LD[[2]]

LD_r<-LD[[3]]
}

