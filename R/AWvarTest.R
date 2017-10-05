# modified on Oct. 21, 2016
#  (1) copy the content of 'iAWvar.c' to 'AWvarTest' and delete 'iAWvar.c'
#   since iAWvar.c is actually AW's equal variance test
#
# modified on Dec. 13, 2015
#  (1) rename 'x' to 'value'
#  (2) rename 'y' to 'group'
#  (3) rename 'scoreTestVarAW.default' to 'AWvarTest'
# created on July 17, 2015
# test for equal variance based on score test of logistic regression

# group - vector of binary values
# value - continuous variable
AWvarTest=function(value,group)
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
  
  # within-group mean centered
  pos1=which(group==1)
  pos0=which(group==0)
  value1=value[pos1]
  value0=value[pos0]
  
  m.value1=mean(value1, na.rm=TRUE)
  m.value0=mean(value0, na.rm=TRUE)
  
  value1.c=(value1-m.value1)^2
  value0.c=(value0-m.value0)^2
  
  z=rep(NA, length(value))
  z[pos1]=value1.c
  z[pos0]=value0.c
  zbar=mean(z, na.rm=TRUE)
  
  groupbar=mean(group, na.rm=TRUE)
  
  U2=sum((group-groupbar)*z, na.rm=TRUE)
  varU2 = groupbar*(1-groupbar)*sum((z-zbar)^2, na.rm=TRUE)
  
  T2 = U2^2/varU2
  pval= 1-pchisq(T2, df=1)
  
  res=list(U2=U2, varU2=varU2, stat=T2, pval=pval, z=z, zbar=zbar) 
  return(res)

}


