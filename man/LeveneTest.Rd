\name{LeveneTest}
\alias{LeveneTest}
\title{
Test for equality of variance based on Levene test
}
\description{
Test for equality of variance based on Levene test.
}
\usage{
LeveneTest(value, group)
}
\arguments{
  \item{value}{
numeric. Measurements to be compared between two groups.
}
  \item{group}{
numeric. Subject's group membership. Must be binary (i.e., taking values 0 or 1).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with 2 elements:
  \item{stat}{test statistic value}
  \item{pval}{pvalue of the score test}
}
\references{
Levene H (1960)
Robust tests for equality of variances.
Contributions to probability and statistics: Essays in honor of Harold Hotelling, 2, 278-292.

Li X, Qiu W, Morrow J, DeMeo DL, Weiss ST, Fu Y, Wang X. (2015) 
A Comparative Study of Tests for Homogeneity of Variances with Application to DNA Methylation Data. 
PLoS ONE 10(12): e0145295. PMID: 26683022

Qiu W, Li X, Morrow J, DeMeo DL, Weiss ST, Wang X, Fu Y. 
New Score Tests for Equality of Variances in the Application of DNA Methylation Data Analysis [Version 2]. 
Insights Genet Genomics.  (2017) 1: 3.2

Li X, Qiu W, Fu Y, Wang X. (2017)
Robust Joint Score Tests in the Application of DNA Methylation Data Analysis.
In submission.

}
\author{
Xuan Li <lixuan0759@mathstat.yorku.ca>,
Weiliang Qiu <stwxq@channing.harvard.edu>,
Yuejiao Fu <yuejiao@mathstat.yorku.ca>,
Xiaogang Wang <stevenw@mathstat.yorku.ca>
}
%\note{
%%  ~~further notes~~
%}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
    # generate simulated data set from t distribution
    set.seed(1234567)
    es.sim = genSimData.tDistr(nCpGs = 100, nCases = 20, nControls = 20,
      df0 = 10, ncp0 = 0, df1 = 6, ncp1 = 2.393, testPara = "var",
      eps = 1.0e-3, applier = lapply) 
    print(es.sim)
    print(exprs(es.sim)[1:2,1:3])

    # do AW score test for the first probe
    dat = exprs(es.sim)
    pDat = pData(es.sim)
    print(pDat[1:2,])

    res = LeveneTest(value = dat[1,], group = pDat$memSubj)
    print(names(res))
    print(res)
}
\keyword{ method }
