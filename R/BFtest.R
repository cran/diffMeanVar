# renamed on Dec. 11, 2015
#  (1) rename 'x' by 'value'
#  (2) rename 'y' by 'group'
#
# renamed on Dec. 3, 2015
#  (1) rename from BFWrapper to BFtest
#
# modified on July 24, 2015
#  (1) write my own function to perform Levene's test since levene.test() is too slow
#  (2) for Brown and Forsythe's test
#
# created on July 23, 2015
# a wrapper function for 
#  Levene's test to test equality of variance
#   for whole genome
#

# test for equal variance
# group - vector of binary values
# value - continuous variable
BFtest=function(value,group)
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

  value1.m=median(value1, na.rm=TRUE)
  value0.m=median(value0, na.rm=TRUE)

  z1=abs(value1 - value1.m)
  z0=abs(value0 - value0.m)

  z=c(z1, z0)
  zbar=mean(z, na.rm=TRUE)
  z1.m=mean(z1, na.rm=TRUE)
  z0.m=mean(z0, na.rm=TRUE)

  n1=length(pos1)
  n0=length(pos0)

  numer=n1*(z1.m-zbar)^2 + n0*(z0.m-zbar)^2
  denom = sum( (z1-z1.m)^2, na.rm=TRUE )
  denom = denom + sum( (z0-z0.m)^2, na.rm=TRUE )

  N=length(value)

  W = (N-2)*numer/ denom
  
  pval = 1 - pf(W, df1=1, df2=N-2)
  
  res=list(stat=W, pval=pval) 

  return(res)
}


