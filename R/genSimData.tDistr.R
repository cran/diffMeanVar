#  revised on July 7, 2016
#  (1) add functions 'dfNcptDistr', 'rootFunc', 'fFunc' to obtain
#      df and ncp for a non-central t distribution given its mean and sd
#      based on the approximation formula shown in
#      https://en.wikipedia.org/wiki/Cubic_function
#
#  (2) add function 'genSimData.tDistr2' to generate simulated data by
#      providing mean and sd (instead of df and ncp) 
#      of non-central t distribution 
#
# revised on July 19, 2015
#  (1) add outlierFlag 
# created on July 13, 2015
#   simulate data set with so that cases and controls
#     have different means and variances


# generate simulated data using the model (Table 1) used in Ahn and Wang (2013)
#  Pacific Symposium on Biocomputing. pp69-79
#
#  T = (Z+mu)/sqrt(V/v) ~ noncentral t with df v and noncentrality parameter mu
#     Z~N(0, 1), V~chisq with df v, Z and V independent
#  
#  E(T) =mu*sqrt(v/2) *Gamma( (v-1)/2 ) / Gamma(v/2), for v>1 
#  var(T) = v*(1+mu^2)/(v-2) - mu^2*v/2*( Gamma( (v-1)/2 ) / Gamma( v/2 ))^2 for v>2
#
# controls are from t10 (mean=0, var=df/(df-2)=10/8=1.25)
# cases are from t.(df=6, ncp=2.393)  (mean=2.754923, var=2.50007)

# ncp=noncentrality parameter mu
# df=degree of freedom
meanVartDistr=function(df, ncp)
{
  if(df<=2)
  {
    stop("df should be > 2!\n")
  }
  if(ncp <0)
  {
    stop("ncp should be >=0!\n")
  }
  mu=ncp
  v=df

  if(abs(mu)<1.0e-6)
  {
    ET=0
    if(v>2)
    {
      varT = v/(v-2)
    } else {
      varT = NA
      cat("if ncp=0, then  df should > 2\n")
    }
  } else {

    if(v>1)
    {
      ET =mu*sqrt(v/2) *gamma( (v-1)/2 ) / gamma(v/2)
    } else {
      ET = NA
      cat("warning! df<=1 so that E(T)=NA\n")
    }
    if(v>2)
    {
      varT = v*(1+mu^2)/(v-2) - mu^2*v/2*( gamma( (v-1)/2 ) / gamma( v/2 ))^2
      if(varT<= 0)
      {
        varT = NA
        cat("warning! varT <=0\n")
      }
    } else {
      varT = NA
      cat("warning! df<=2 so that var(T)=NA\n")
    }
  }
 
  res=c(ET, varT)
  names(res)=c("mean", "variance")
  return(res)
}

##########################################
# An approximate formulas
#   theta = mu / [1 - 3/(4*nu-1)]
#   sigma2 = nu/(nu-2)*(1+mu^2)-mu^2/[1-3/(4*nu-1)]^2
#  where theta is the mean, sigma2 is the variance
#  nu is the degree of freedom, mu is the non-centrality parameter
# https://en.wikipedia.org/wiki/Noncentral_t-distribution
#
# We can then get
#   (a) nu = (4*theta-mu)/[4*(theta-mu)]
#   (b) mu^3-4*theta*mu^2+mu*[7*(sigma^2+theta^2)+1] -
#       4*theta*(sigma^2+theta^2+1) = 0

# The critical point for the function
#   f(x)=a*x^3+b*x^2+c*x+d
#  is
#  x_{critical}=[-b \pm \sqrt{b^2-3*a*c}]/(3*a)
#
# critical points are values of x where the slope of the function
#   is zero
#
# the inflection point is
#  x_{inflection} = -b/(3*a)
#
# if b^2-3*a*c<0, then f(x) is monotonic.
#  Hence, cubic equation f(x)=0 has one and only one unique solution.
#
# for our case, a=1, b=-4*theta, c=7*(sigma^2+theta^2)+1
#  d = -4*theta*(sigma^2+theta^2+1)
# then
# b^2-3*a*c=-21*sigma^2-5*theta^2-3<0

# theta - mean of noncentral t distribution
# sigma - standard deviation of non-central t distribution 
# obtain root for cubic equation
rootFunc=function(theta, sigma)
{

  a=1
  b=-4*theta
  myc=7*(sigma^2+theta^2)+1
  d= -4*theta*(sigma^2+theta^2+1)

  delta0=b^2-3*a*myc
  delta02=-21*sigma^2-5*theta^2-3
  #cat("delta0=", delta0, "delta02=", delta02, "\n")

  delta1=2*b^3-9*a*b*myc+27*a^2*d
  delta1.2=theta*(16*theta^2+144*sigma^2-72)

  #cat("delta1=", delta1, ", delta1.2=", delta1.2, "\n")

  part1=sqrt(delta1.2^2-4*delta02^3)
  part2.1=delta1.2+part1
  part2.2=delta1.2-part1

  CaptialC1=( part2.1/2 )^{1/3}
  #CaptialC2=( part2.2/2 )^{1/3}
  #cat("CaptialC1=", CaptialC1, ", CaptialC2=", CaptialC2, "\n")

  x = b+CaptialC1+delta02/CaptialC1
  x = - x / (3*a)

  #cat("x=", x, "\n")
  return(x)
}


# theta - mean of non-central t distribution
# sigma - standard deviation of non-central t distribution

# given mean and sd of non-central t distribution, obtain
#  df and ncp
dfNcptDistr=function(theta, sigma)
{
  if(sigma<=1)
  {
    stop("sigma should be >1!\n")
  }
  # if mean = 0
  if(abs(theta)<1.0e-6)
  {
    mu.opt=0
    if(sigma^2>1)
    {
      nu.opt=2*sigma^2/(sigma^2-1)
    } else {
      nu.opt=NA
      cat("warning! if mean=0, then variance sigma^2 should be > 1!\n")
    }
  } else {
    mu.opt=rootFunc(theta=theta, sigma=sigma)
    nu.opt=(4*theta-mu.opt)/(4*(theta-mu.opt))
  }
  
  res=c(df=nu.opt, ncp=mu.opt)
  return(res)
}

genSimData.tDistr.default=function(nCases, nControls,
  df0=10, ncp0=0, df1=6, ncp1=2.393) 
{
  controls=rt(n=nControls, df=df0, ncp=ncp0)
  cases=rt(n=nCases, df=df1, ncp=ncp1)

  x=c(cases, controls)
  y=c(rep(1, nCases), rep(0, nControls))
  frame=data.frame(x=x, y=y)
  invisible(frame)
}

genSimData.tDistr=function(nCpGs,nCases, nControls,
  df0=10, ncp0=0, df1=6, ncp1=2.393, testPara="var", 
  outlierFlag = FALSE,
  eps=1.0e-3, applier=lapply) 
{
  y=c(rep(1, nCases), rep(0, nControls))
  ttLst=applier(1:nCpGs, function(i) {
    resi=genSimData.tDistr.default(nCases=nCases, nControls=nControls,
      df0=df0, ncp0=ncp0, df1=df1, ncp1=ncp1) 
    return(resi$x)
  })
  mat=t(sapply(ttLst, function(x) { x} ))

  if(outlierFlag)
  {
    # add outlier
    datvec=c(mat)
    M.max=max(datvec, na.rm=TRUE)
  
    # randomly select one case to add outliers
    pos.outCase=sample(x=1:nCases, size=1, replace=FALSE)
    mat[, pos.outCase]=M.max
  }
 
  mean.var0=meanVartDistr(df=df0, ncp=ncp0)
  mean.var1=meanVartDistr(df=df1, ncp=ncp1)

  m0=mean.var0[1]
  v0=mean.var0[2]

  m1=mean.var1[1]
  v1=mean.var1[2]

  if(testPara=="mean")
  {
     if(abs(m0-m1)<eps) 
     {
       memGenes=rep(0, nCpGs)
     } else {
       memGenes=rep(1, nCpGs)
     }
  } else if(testPara=="var") {
     if(abs(v0-v1)<eps) 
     {
       memGenes=rep(0, nCpGs)
     } else {
       memGenes=rep(1, nCpGs)
     }
  } else { # test for both mean difference and variance difference
     if(abs(m0-m1)<eps && abs(v0-v1)<eps) 
     {
       memGenes=rep(0, nCpGs)
     } else {
       memGenes=rep(1, nCpGs)
     }
  }
   
  nSubj=nCases+nControls
  subjID=paste("subj", 1:nSubj, sep="")
  probeID=paste("probe", 1:nCpGs, sep="")
  genes=paste("gene", 1:nCpGs, sep="")
  chr=rep(1, nCpGs)

  pDat=data.frame(
    arrayID=subjID,    
    memSubj=y
  )
  rownames(pDat)=subjID

  fDat=data.frame(probe=probeID, gene=genes, chr=chr, memGenes=memGenes)
  rownames(fDat)=probeID

  rownames(mat)=probeID
  colnames(mat)=subjID

  # create ExpressionSet object
  aa<-new("AnnotatedDataFrame", data=pDat)
  bb = as(mat, "matrix")

  es<-new("ExpressionSet",
      exprs=bb,
      phenoData = aa,
      annotation = "")

  Biobase::fData(es)=fDat

  invisible(es)
}


genSimData.tDistr2=function(nCpGs,nCases, nControls,
  mu0=0, sigma0=sqrt(1.25), mu1=2.75, sigma1=sqrt(2.5), 
  testPara="var", 
  outlierFlag = FALSE,
  eps=1.0e-3, applier=lapply) 
{
  res0=dfNcptDistr(theta=mu0, sigma=sigma0)
  res1=dfNcptDistr(theta=mu1, sigma=sigma1)

  df0=res0[1]
  ncp0=res0[2]

  df1=res1[1]
  ncp1=res1[2]

  y=c(rep(1, nCases), rep(0, nControls))
  ttLst=applier(1:nCpGs, function(i) {
    resi=genSimData.tDistr.default(nCases=nCases, nControls=nControls,
      df0=df0, ncp0=ncp0, df1=df1, ncp1=ncp1) 
    return(resi$x)
  })
  mat=t(sapply(ttLst, function(x) { x} ))

  if(outlierFlag)
  {
    # add outlier
    datvec=c(mat)
    M.max=max(datvec, na.rm=TRUE)
  
    # randomly select one case to add outliers
    pos.outCase=sample(x=1:nCases, size=1, replace=FALSE)
    mat[, pos.outCase]=M.max
  }
 
  mean.var0=meanVartDistr(df=df0, ncp=ncp0)
  mean.var1=meanVartDistr(df=df1, ncp=ncp1)

  m0=mean.var0[1]
  v0=mean.var0[2]

  m1=mean.var1[1]
  v1=mean.var1[2]

  if(testPara=="mean")
  {
     if(abs(m0-m1)<eps) 
     {
       memGenes=rep(0, nCpGs)
     } else {
       memGenes=rep(1, nCpGs)
     }
  } else if(testPara=="var") {
     if(abs(v0-v1)<eps) 
     {
       memGenes=rep(0, nCpGs)
     } else {
       memGenes=rep(1, nCpGs)
     }
  } else { # test for both mean difference and variance difference
     if(abs(m0-m1)<eps && abs(v0-v1)<eps) 
     {
       memGenes=rep(0, nCpGs)
     } else {
       memGenes=rep(1, nCpGs)
     }
  }
   
  nSubj=nCases+nControls
  subjID=paste("subj", 1:nSubj, sep="")
  probeID=paste("probe", 1:nCpGs, sep="")
  genes=paste("gene", 1:nCpGs, sep="")
  chr=rep(1, nCpGs)

  pDat=data.frame(
    arrayID=subjID,    
    memSubj=y
  )
  rownames(pDat)=subjID

  fDat=data.frame(probe=probeID, gene=genes, chr=chr, memGenes=memGenes)
  rownames(fDat)=probeID

  rownames(mat)=probeID
  colnames(mat)=subjID

  # create ExpressionSet object
  aa<-new("AnnotatedDataFrame", data=pDat)
  bb = as(mat, "matrix")

  es<-new("ExpressionSet",
      exprs=bb,
      phenoData = aa,
      annotation = "")

  Biobase::fData(es)=fDat

  invisible(es)
}




