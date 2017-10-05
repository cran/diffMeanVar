# modified on Dec. 13, 2015
#  (1) rename 'y' to 'group'
#  (2) rename 'x' to 'value'
#  (3) rename 'FtestWrapper' to 'FTest'
#
# created on July 23, 2015
# a wrapper function for 
#  F test to test equality of variance
#   for whole genome
#

# test for equal variance
# group - vector of binary values
# value - continuous variable
FTest=function(value,group)
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

  res.L=stats::var.test(x=value[which(group==0)], y=value[which(group==1)])

  res=list(stat=res.L$statistic, pval=res.L$p.value) 
  return(res)
}


