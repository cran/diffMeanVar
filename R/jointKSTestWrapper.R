# created on Dec. 13, 2015
#  (1) ks test for equal distribution
#
# created on Dec. 13, 2015
#  (1) rename 'LRTequalMeanVarWrapper' to 'jointLRTTestsWrapper'
#  (2) rename 'LRTequalMeanVar' to 'jointLRTTests'
#  (3) rename 'D' in the output of 'jointLRTTests' to 'stat'
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
jointKSTest=function(value,group)
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

  res.ks=stats::ks.test(x=value[which(group==0)], y=value[which(group==1)])
  res=list(stat=res.ks$statistic, pval=res.ks$p.value) 

  return(res)

}

# esFlag = "es" (for ExpressionSet object) or "meth" (for methylumi object)
jointKSTestWrapper=function(
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
    ans=try(resi <- jointKSTest(value=xi,group=y))
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

