\name{permutation}
\title{permutation for case-control data}
\description{
This function is to permute cases and controls randomly.}

\usage{ permuC(phenotype,genotype) }
\arguments{
\item{phenotype}{ numeric vector with phenotype status: 0=controls, 1=cases. No missing data allowed }
\item{genotype}{ numeric matrix or data frame with genotype data coded as 0, 1, 2. No missing data is allowed.}
}

\value{The matrix in which the first column is the phenotype and others are the corresponding genotype after permutation.}

\references{Chan Wang, Leiming Sun, Liming Li and Yue-Qing Hu, Retrospective Score Test for Case-Control Association Study}
\author{Chan Wang and Yue-Qing Hu }

