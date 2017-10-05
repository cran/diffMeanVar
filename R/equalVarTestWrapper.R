# wrapper function to test for equality of variance

# esFlag = "es" (for ExpressionSet object) or "meth" (for methylumi object)
# varTestFunc - user defined function with at least 2 arguments: value and group
#  possible varTestFunc provided by 'varianceDiff' package are:
#   'iAWvar.c', 'iAWvar.BF', 'iAWvar.FK', 'iAWvar.Levene', 'iAWvar.TrimMean'
#   'LRTequalVar',  'AWvarTest',
#   'BFTest', 'FTest', 'BartlettTest'
#   'LeveneTest', 'FKTest', 'TrimMeanLeveneTest', 'flignerTest'
equalVarTestWrapper=function(
  es,
  grpVar = "group", 
  varTestFunc = AWvarTest,
  esFlag = "es", 
  pvalAdjMethod = "fdr",
  alpha = 0.05,
  nTop = 20,
  probeID.var = "ProbeID",
  gene.var = "Symbol",
  chr.var = "Chromosome",
  applier = lapply,
  verbose = FALSE)
{

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
    ans=try(resi <- varTestFunc(value=xi,group=y))
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
  frame=data.frame(probe=probeID, gene=genes, chr=chr,
    pos=1:nr, stat=mat[,1], pval=mat[,2])
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

