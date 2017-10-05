\name{genSimData.tDistr}
\alias{genSimData.tDistr}
\title{
Generating simulated data set from t distributions
}
\description{
Generating simulated data set from t distributions.
}
\usage{
genSimData.tDistr(
  nCpGs, 
  nCases, 
  nControls,
  df0 = 10, 
  ncp0 = 0, 
  df1 = 6, 
  ncp1 = 2.393, 
  testPara = "var", 
  outlierFlag = FALSE,
  eps = 0.001, 
  applier = lapply)
}
\arguments{
  \item{nCpGs}{
integer. Number of genes.
}
  \item{nCases}{
integer. Number of cases.
}
  \item{nControls}{
integer. Number of controls.
}
  \item{df0}{
integer. Degree of freedom for group 0 (controls).
}
  \item{ncp0}{
numeric. Non-centrality parameter for group 0 (controls).
}
  \item{df1}{
integer. Degree of freedom for group 1 (cases).
}
  \item{ncp1}{
numeric. Non-centrality parameter for group 1 (cases).
}
  \item{testPara}{
character string. Indicating if the test is for testing equal mean 
(\code{testPara="mean"}), for testing equal variance
(\code{testPara="var"}), or for testing both equal mean
and equal variance (\code{testPara="both"}).
}
  \item{outlierFlag}{
logical. Indicating if outliers would be generated. If \code{outlierFlag=TRUE},
then we followed Phipson and Oshlack's (2014) simulation studies to 
generate one outlier for each CpG site by replacing the DNA methylation 
level of one diseased subject by the maximum of the DNA methylation levels 
of all CpG sites.
}
  \item{eps}{
numeric. If \eqn{|mean_0-mean_1|<eps}{|mean0-mean1|<eps} then we regard \eqn{mean_0=mean_1}{mean0=mean1}.
Similarly, 
if \eqn{|var_0-var_1|<eps}{|var0-var1|<eps} then we regard \eqn{var_0=var_1}{var0=var1}.
\eqn{mean_0}{mean0} and \eqn{var_0}{var0} are the mean and variance of
the t distribution for controls.
\eqn{mean_1}{mean1} and \eqn{var_1}{var1} are the mean and variance of
the t distribution for cases.
}
  \item{applier}{
function name to do \code{apply} operation.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An ExpressionSet object. The phenotype data of the ExpressionSet object
contains 2 columns: \code{arrayID} (array id) and {memSubj} (subject
membership, i.e., case (\code{memSubj=1}) or control (\code{memSubj=0})).
The feature data of the ExpressionSet object contains 4 elements:
\code{probe} (probe id), \code{gene} (psuedo gene symbol), \code{chr}
(psuedo chromosome number), and \code{memGenes}
(indicating if a gene is differentially expressed (when \code{testPara="mean"})
or indicating if a gene is differentially variable (when \code{testPara="var"}) ).

}
\references{
  Ahn S. and Wang T. (2013)
  A Powerful Statistical Method for Indentifying Differentially Methylated Markers in Complex Diseases.
  Pacific Symposium on Biocomputing. 69-79.


  Phipson B, Oshlack A. 
  DiffVar: A new method for detecting differential variability with application to methylation in cancer and aging. 
  \emph{Genome Biol} 2014; 15:465

  Li X, Qiu W, Morrow J, DeMeo DL, Weiss ST, Fu Y, Wang X. (2015) 
  A Comparative Study of Tests for Homogeneity of Variances with Application to DNA Methylation Data. 
  PLoS ONE 10(12): e0145295. PMID: 26683022
  
  Qiu W, Li X, Morrow J, DeMeo DL, Weiss ST, Wang X, Fu Y. 
  New Score Tests for Equality of Variances in the Application of DNA Methylation Data Analysis [Version 2]. 
  Insights Genet Genomics.  (2017) 1: 3.2
  

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
  %\dontrun{
    # generate simulated data set from t distribution
    set.seed(1234567)
    es.sim = genSimData.tDistr(nCpGs = 100, nCases = 20, nControls = 20,
      df0 = 10, ncp0 = 0, df1 = 6, ncp1 = 2.393, testPara = "var",
      outlierFlag = FALSE, 
      eps = 1.0e-3, applier = lapply) 
    print(es.sim)
  %}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
