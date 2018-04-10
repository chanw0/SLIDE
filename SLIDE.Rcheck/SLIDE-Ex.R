pkgname <- "SLIDE"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('SLIDE')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("LD-EM")
### * LD-EM

flush(stderr()); flush(stdout())

### Name: LD-EM
### Title: Estimating LD using EM algorithm

### ** Examples

genotype<-matrix(sample(c(0,1,2),5000,replace=TRUE),500,10) ### 500 individuals and 10 SNPs

LD<-ld.em(genotype)

LD_orignal<-LD[[1]]

LD_D<-LD[[2]]

LD_r<-LD[[3]]



cleanEx()
nameEx("SLIDE")
### * SLIDE

flush(stderr()); flush(stdout())

### Name: SLIDE
### Title: Retrospective Score Test for Case-Control Association Study

### ** Examples


library(MASS)

set.seed(1234)

genotype<-matrix(sample(c(0,1,2),5000,replace=TRUE),500,10) ### 500 individuals and 10 SNPs

phenotype<-c(rep(1,200),rep(0,300)) ### 200 cases and 300 controls

Statistics<-SILDE(cbind(phenotype,genotype),method=c("Cov"),LD=NULL,num_per=200)

SLIDE<-Statistics[[1]]

## 21.2053, 0.0200, 0.0197




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
