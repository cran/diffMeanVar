# modified on July 15, 2016
#  (1) revise the definition of zi
#      zi1=|xi1-trimedmean(xi1)|
#      zi0=|xi0-trimedmean(xi0)|
#
# modified on Dec. 13, 2015
#  (1) rename 'y' to 'group'
#  (2) rename 'x' to 'value'
#  (3) rename scoreTestVarTrimMean.default' to 'iAWvar.TrimMean'
#
# created on July 19, 2015
# scor test of logistic regresssion based on Trimmed-mean based 
#  Levene's test to test equality of variance
#


# test for equal variance
# group - vector of binary values
# value - continuous variable
# trim.alpha=0.25 is default value for R function 'mean', it is
#  the fraction (0 to 0.5) of observations to be trimmed from each 
#  end of x before the mean is computed. Values of trim.alpha outside 
#  that range are taken as the nearest endpoint.
iAWvar.TrimMean=function(value,group, trim.alpha = 0.25)
{
  u.group=sort(unique(group))
  if(length(u.group)!=2)
  {
    stop("group must take 2 and only 2 values\n")
  }
  if(!identical(u.group, c(0, 1)))
  {
    stop("group must only take values 0 or 1\n")
  }
  if(length(value) != length(group))
  {
    stop("value must have the same length as group\n")
  }

  pos1=which(group==1)
  pos0=which(group==0)

  value1=value[pos1]
  value0=value[pos0]

  # get trimmed mean
  m.value1=mean(value1, na.rm=TRUE, trim = trim.alpha)
  m.value0=mean(value0, na.rm=TRUE, trim = trim.alpha)

  # trimmed mean centering
  value1.2=abs(value1-m.value1)
  value0.2=abs(value0-m.value0)

  z=rep(NA, length(value))
  z[pos1]=value1.2
  z[pos0]=value0.2

  ybar=mean(group, na.rm=TRUE)

  U2=sum((group-ybar)*z, na.rm=TRUE)

  zbar=mean(z, na.rm=TRUE)

  varU2 = ybar*(1-ybar)*sum((z-zbar)^2, na.rm=TRUE)

  T2 = U2^2/varU2
  pval= 1-pchisq(T2, df=1)

  res=list(U2=U2, varU2=varU2, stat=T2, pval=pval, z=z, zbar=zbar) 
  return(res)
}


