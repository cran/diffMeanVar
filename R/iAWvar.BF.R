# revised on July 15, 2016
#  (1) use new definition: 
#      zi1=|xi1-median(xi1)|
#      zi0=|xi0-median(xi0)|
#
# revised on Dec. 12, 2015
#  (1) replace 'y' by 'group'
#  (2) replace 'x' by 'value'
# renamed on Dec. 3, 2015
#  rename scoreTestVarBF.default to iAWvar.BF
#
# created on July 17, 2015
# scor test of logistic regresssion based on Brown and Forsythe's test
# to test equality of variance
#


# test for equal variance
# group - vector of binary values
# value - continuous variable
iAWvar.BF=function(value,group)
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

  # get median
  m.value1=median(value1, na.rm=TRUE)
  m.value0=median(value0, na.rm=TRUE)

  # median centering
  value1.2=abs(value1-m.value1)
  value0.2=abs(value0-m.value0)

  z=rep(NA, length(value))
  z[pos1]=value1.2
  z[pos0]=value0.2

  groupbar=mean(group, na.rm=TRUE)

  U2=sum((group-groupbar)*z, na.rm=TRUE)

  zbar=mean(z, na.rm=TRUE)

  varU2 = groupbar*(1-groupbar)*sum((z-zbar)^2, na.rm=TRUE)

  T2 = U2^2/varU2
  pval= 1-pchisq(T2, df=1)

  res=list(U2=U2, varU2=varU2, stat=T2, pval=pval, z=z, zbar=zbar) 
  return(res)
}


