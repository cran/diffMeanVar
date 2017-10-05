# revised on Dec. 12, 2015
#  (1) rename 'y' to 'group'
#  (2) rename 'x' to 'value'
#  (3) rename 'scoreTestMean.default' to 'scoreTestMean'
#
# created on July 17, 2015
# test for equal mean based on score test of logistic regression

# group - vector of binary values
# value - continuous variable
scoreTestMean=function(value,group)
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

  # obtain U1
  groupbar=mean(group, na.rm=TRUE)
  valuebar=mean(value, na.rm=TRUE)
  U1=sum((group-groupbar)*value, na.rm=TRUE)

  varU1 = groupbar*(1-groupbar)*sum((value-valuebar)^2, na.rm=TRUE)

  T1 = U1^2/varU1
  pval= 1-pchisq(T1, df=1)

  res=list(U1=U1, varU1=varU1, stat=T1, pval=pval, x=value, xbar=valuebar) 
  return(res)
}


