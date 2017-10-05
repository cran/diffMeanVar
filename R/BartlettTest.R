# modified on Dec. 12, 2015
#  (1) rename 'y' to 'group'
#  (2) rename 'x' to 'value'
#  (3) rename 'BartlettWrapper' to 'BartlettTest'
#
# created on July 23, 2015
# a wrapper function for 
#  Bartlett's test to test equality of variance
#   for whole genome
#

# test for equal variance
# group - vector of binary values
# value - continuous variable
BartlettTest=function(value,group)
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

  ttframe=data.frame(x=value, memSubj=group)
  res.B=stats::bartlett.test(x~memSubj, data=ttframe)

  res=list(stat=res.B$statistic, pval=res.B$p.value) 
  return(res)
}


