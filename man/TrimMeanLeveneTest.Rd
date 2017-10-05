\name{TrimMeanLeveneTest}
\alias{TrimMeanLeveneTest}
\title{
Test for equality of variance based on trimmed-mean based Levene test
}
\description{
Test for equality of variance based on trimmed-mean based Levene test.
}
\usage{
TrimMeanLeveneTest(value, group, trim.alpha = 0.25)
}
\arguments{
  \item{value}{
numeric. Measurements to be compared between two groups.
}
  \item{group}{
numeric. Subject's group membership. Must be binary (i.e., taking values 0 or 1).
}
  \item{trim.alpha}{
numeric. Indicating proportion of data points to be trimmed from both ends before calculating sample mean.
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
Brown MB and Forsythe AB (1974)
Robust Tests for Equality of Variances.
Journal of the American Statistical Association, 69, 364-367.


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

    res = TrimMeanLeveneTest(value = dat[1,], group = pDat$memSubj)
    print(names(res))
    print(res)
}
\keyword{ method }
