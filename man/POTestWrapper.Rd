\name{POTestWrapper}
\alias{POTestWrapper}
\title{
Wrapper function to test for equality of variance by using Phipson and Oshlack's (2014) methods
}
\description{
Wrapper function to test for equality of variance by using Phipson and Oshlack's (2014) methods
}
\usage{
POTestWrapper(
  es, 
  grpVar = "group", 
  type = "AD",
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
  \item{type}{
character string indicating if \code{AD} method or \code{SQ} method
would be used.
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
determine if a probe is diferentially variable.
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
logical. indicating if intermediate results should be output to screen.
}
}
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
Phipson B, Oshlack A.
DiffVar: a new method for detecting differential variability with application to methylation in cancer and aging.
\emph{Genome Biology} 2014, 15:465.
}
\author{
Xuan Li <lixuan0759@mathstat.yorku.ca>,
Weiliang Qiu <stwxq@channing.harvard.edu>,
Yuejiao Fu <yuejiao@mathstat.yorku.ca>,
Xiaogang Wang <stevenw@mathstat.yorku.ca>
}

\examples{
    # generate simulated data set from t distribution
    set.seed(1234567)
    es.sim = genSimData.tDistr(nCpGs = 100, nCases = 20, nControls = 20,
      df0 = 10, ncp0 = 0, df1 = 6, ncp1 = 2.393, testPara = "var",
      eps = 1.0e-3, applier = lapply) 
    print(es.sim)
    
    res.POTestWrapper = POTestWrapper(
      es = es.sim,
      grpVar = "memSubj",
      type = "SQ",
      esFlag = "es", 
      pvalAdjMethod = "fdr",
      alpha = 0.05,
      nTop = 20,
      probeID.var = "probe",
      gene.var = "gene",
      chr.var = "chr",
      applier=lapply,
      verbose=TRUE)

}
\keyword{ method }
