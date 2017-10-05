\name{AWvarTest}
\alias{AWvarTest}
\title{
Test for equality of variance based on Ahn and Wang's (2013) score test
}
\description{
Test for equality of variance based on Ahn and Wang's (2013) score test.
}
\usage{
AWvarTest(value, group)
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
  \item{pval}{p-value of the score test}
  \item{z}{squared deviation of \code{value} from mean \code{value}}
  \item{zbar}{sample average of \code{z}}
}
\references{
Ahn S. and Wang T. (2013)
A Powerful Statistical Method for Indentifying Differentially Methylated Markers in Complex Diseases.
Pacific Symposium on Biocomputing. 69-79.
 
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

    res = AWvarTest(value = dat[1,], group = pDat$memSubj)
    print(names(res))
    print(res)
}
\keyword{ method }
