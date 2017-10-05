\name{equalMeanTestWrapper}
\alias{equalMeanTestWrapper}
\title{
Wrapper function to test for equality of mean by using the score test of logistic regression for whole genome
}
\description{
Wrapper function to test for equality of mean by using the score test of logistic regression for whole genome.
}
\usage{
equalMeanTestWrapper(
  es, 
  grpVar = "group", 
  meanTestFunc = scoreTestMean,
  esFlag = "es", 
  pvalAdjMethod = "fdr", 
  alpha = 0.05, 
  nTop = 20, 
  probeID.var = "ProbeID", 
  gene.var = "Symbol", 
  chr.var = "Chromosome", 
  applier = lapply, 
  verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{es}{
An ExpressionSet object storing gene expression/DNA methylation data,
phenotype data, and feature annotation.
}
  \item{grpVar}{
character string. The name of the phenotype variable indicating arrays' group
membership. 0 means control and 1 means case.
}
  \item{meanTestFunc}{
A user-defined function to test for equality of mean. 
Available functions
provided by \code{diffMeanVar} package include:
\code{scoreTestMean}. 
The function
should contain 2 and only 2 inputs: \code{value} and \code{group},
where \code{group} is a binary vector indicating if a subject
is a case (\code{group=1}) or control (\code{group=0}) and \code{value}
is a continuous-type variable. The goal of the function 
is to test if the mean of \code{value} for cases
is the same as that for controls. The function should return
a list object containing at least 2 elements: \code{stat} and \code{pval}. 
For each probe stored in \code{es}, the function will be applied
to test if the probe is different in means between cases and
controls.
}
  \item{esFlag}{
character string. Indicating if \code{es} is an \code{ExpressionSet} object 
or \code{MethylSet} object. The program will use \code{exprs} function
to extract gene expression data or use \code{betas} function
to extract methylation data.
}
  \item{pvalAdjMethod}{
character string. Indicating which p-value adjustment will be used
to control for multiple testing.
}
  \item{alpha}{
numeric.  Cutoff for p-value or adjusted p-value to
determine if a probe is diferentially mean.
}
  \item{nTop}{
integer.  Specifying the number of top probes to be displayed
if \code{verbose=TRUE}
}
  \item{probeID.var}{
character string. Feature annotation variable indicating probe ID.
}
  \item{gene.var}{
character string. Feature annotation variable indicating gene symbol.
}
  \item{chr.var}{
character string. Feature annotation variable indicating chromosome number.
}
  \item{applier}{
function name to do \code{apply} operation.
}
  \item{verbose}{
logical. Indicating if intermediate results should be output to screen.
}
}
%\details{
%}
\value{
A list of 2 elements. The first element \code{frame} is unsorted data frame;
the second element \code{frame.s} is a sorted data frame object storing the analysis results
and containing the following columns:
\code{probe} (probe id), \code{stat} (test statistic), 
\code{pval} (raw p-value), \code{p.adj} (adjusted p-value),
\code{gene} (gene symbol), \code{chr} (chromosome number), and 
\code{pos} (position of a probe in unsorted data frame).

The data frame  is sorted based on the descending order of
the absolute value of the test statistic.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
  %\dontrun{
    # generate simulated data set from t distribution
    set.seed(1234567)
    es.sim = genSimData.tDistr(nCpGs = 100, nCases = 20, nControls = 20,
      df0 = 10, ncp0 = 0, df1 = 6, ncp1 = 2.393, testPara = "mean",
      eps = 1.0e-3, applier = lapply) 
    print(es.sim)
    
    res.equalMeanTestWrapper = equalMeanTestWrapper(
      es = es.sim,
      grpVar = "memSubj",
      meanTestFunc = scoreTestMean,
      esFlag = "es", 
      pvalAdjMethod = "fdr",
      alpha = 0.05,
      nTop = 20,
      probeID.var = "probe",
      gene.var = "gene",
      chr.var = "chr",
      applier=lapply,
      verbose=TRUE)

  %}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
