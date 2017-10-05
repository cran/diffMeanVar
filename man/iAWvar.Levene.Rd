\name{iAWvar.Levene}
\alias{iAWvar.Levene}
\title{
Test for equality of variance based on improved Ahn and Wang's (2013) score test
}
\description{
Test for equality of variance based on improved Ahn and Wang's (2013) score test.
}
\usage{
iAWvar.Levene(value, group)
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
A list with 6 elements:
  \item{U2}{score statistic}
  \item{varU2}{estimated variance of the score statistic}
  \item{T2}{score statistic \eqn{U_2^2/var(U_2)}{U2^2/varU2}}
  \item{pval}{pvalue of the score test}
  \item{z}{absolute deviation of \code{value} from mean \code{value}}
  \item{zbar}{sample average of \code{z}}
}
\references{
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
% See Also as \code{\link{help}}
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

    res = iAWvar.Levene(value = dat[1,], group = pDat$memSubj)
    print(names(res))
    print(res)
}
\keyword{ method }
