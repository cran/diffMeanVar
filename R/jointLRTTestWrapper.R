# created on Dec. 13, 2015
#  (1) rename 'LRTequalMeanVarWrapper' to 'jointLRTTestWrapper'
#  (2) rename 'LRTequalMeanVar' to 'jointLRTTest'
#  (3) rename 'D' in the output of 'jointLRTTest' to 'stat'
#
# modified on Dec. 13, 2015
#  (1) rename 'y' to 'group'
#  (2) rename 'x' to 'value'
#  (3) rename 'LRTequalMeanVar' to 'LRTequalMeanVarWrapper'
#  (4) rename 'LRTequalMeanVar.default' to 'LRTequalMeanVar'
#
# created on July 17, 2015
#   jointly test for equal mean or equal variance based on likelihood ratio test

# group - vector of binary values
# value - continuous variable
jointLRTTest=function(value,group)
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

  value1sum=sum(value1, na.rm=TRUE)
  value1.sq.sum=sum(value1^2, na.rm=TRUE)
  n1=length(value1)
  svalue1.sq=(value1.sq.sum-value1sum^2/n1)/n1

  value0sum=sum(value0, na.rm=TRUE)
  value0.sq.sum=sum(value0^2, na.rm=TRUE)
  n0=length(value0)
  svalue0.sq=(value0.sq.sum-value0sum^2/n0)/n0

  n=n1+n0
  p1=n1/n
  p0=n0/n

  m.overall=mean(value, na.rm=TRUE)
  v1.sq=sum((value1-m.overall)^2,na.rm=TRUE)/n1
  v0.sq=sum((value0-m.overall)^2,na.rm=TRUE)/n0

  D=n*log(p1*v1.sq+p0*v0.sq)-n1*log(svalue1.sq)-n0*log(svalue0.sq)

  pval= 1-pchisq(D, df=2)

  res=list(stat=D, pval=pval) 
  return(res)
}

# esFlag = "es" (for ExpressionSet object) or "meth" (for methylumi object)
jointLRTTestWrapper=function(
  es,
  grpVar="group", 
  esFlag="es", 
  pvalAdjMethod="fdr",
  alpha = 0.05,
  nTop=20,
  probeID.var = "ProbeID",
  gene.var = "Symbol",
  chr.var = "Chromosome",
  applier=lapply,
  verbose=FALSE)
{
  # check if  'grpVar' is in 'es' or not

  pDat=Biobase::pData(es)
  fDat=Biobase::fData(es)

  # check if  'grpVar' is in 'es' or not
  y=pDat[, c(grpVar)]

  # number of probes
  nr=nrow(fDat)
  if(!is.null(probeID.var))
  {
    probeID=fDat[, c(probeID.var)]
  } else {
    probeID=paste("pseudoprobe", 1:nr, sep="")
  }

  if(!is.null(gene.var))
  {
    genes=fDat[, c(gene.var)]
  } else {
    genes=paste("gene", 1:nr, sep="")
  }

  if(!is.null(chr.var))
  {
    chr=fDat[, c(chr.var)]
  } else {
    chr=paste("pseudochr", 1:nr, sep="")
  }

  if(esFlag=="es")
  {
    dat=exprs(es)
  } else {
    dat=betas(es)
  }

  ttLst=applier(1:nr, function(i) {
    xi=dat[i,]
    ans=try(resi <- jointLRTTest(value=xi,group=y))
    aaa<-attr(ans,which="class")
    if(!is.null(aaa))
    {
      resi2=rep(NA,2)
    } else {
      resi2=c(resi$stat, resi$pval)   
    }
    names(resi2)=c("stat", "pval")
    return(resi2)
  })
  mat=t(sapply(ttLst, function(x) { return(x) }))
  frame=data.frame(probe=probeID, gene=genes, chr=chr, pos=1:nr, 
    stat=mat[,1], pval=mat[,2])
  frame$p.adj=p.adjust(frame$pval, method=pvalAdjMethod)

  frame.s=frame[order(abs(frame$stat), decreasing = TRUE),]

  if(verbose)
  {
    cat("\n**** output results>>>>\n")
    cat("\n top", nTop, " probes>>\n")
    print(frame.s[1:nTop,])
    cat("\nNumber of probes>>", nr, "\n")
    cat("\nNumber of arrays>>", nrow(pDat), "\n")
    cat("\ngroup variable>>", grpVar, "\n")
    cat("\nfrequencies of groups>>\n")
    print(table(y, useNA="ifany"))
    cat("\nNumber of probes with raw p-value <", alpha, "  >>>", 
      sum(frame.s$pval<alpha, na.rm=TRUE), "\n")
    cat("\nNumber of probes with adjusted p-value <", alpha, " >>>", 
      sum(frame.s$p.adj<alpha, na.rm=TRUE), "\n")
  }

  res=list(frame=frame, frame.s=frame.s)
  invisible(res)
}

