# created on Sept. 9, 2015
#  (1) wrapper function for Phipson and Oshlack's (2014) paper 

# esFlag = "es" (for ExpressionSet object) or "meth" (for methylumi object)
POTestWrapper=function(
  es,
  grpVar = "group", 
  type = "AD",
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
  grp=factor(pDat[, c(grpVar)])

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

  mydesign<-model.matrix(~grp)

  # Fit linear model for differential variability
  # absolute residuals
  if(type=="AD")
  {
    #vfit <- missMethyl::varFit(data=dat, design=mydesign, coef=c(1,2), type="AD")
    vfit <- varFit(data=dat, design=mydesign, coef=c(1,2), type="AD")
  } else { # type=="SQ"
    #vfit <- missMethyl::varFit(data=dat, design=mydesign, coef=c(1,2), type="SQ")
    vfit <- varFit(data=dat, design=mydesign, coef=c(1,2), type="SQ")
  }

  # test statistic
  stat.PO=vfit$t[,2]

  # pvalue
  pval.PO=vfit$p.value[,2]

  frame=data.frame(probe=probeID, gene=genes, chr=chr,
    pos=1:nrow(dat), stat=stat.PO, pval=pval.PO)
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
    print(table(grp, useNA="ifany"))
    cat("\nNumber of probes with raw p-value <", alpha, "  >>>", 
      sum(frame.s$pval<alpha, na.rm=TRUE), "\n")
    cat("\nNumber of probes with adjusted p-value <", alpha, " >>>", 
      sum(frame.s$p.adj<alpha, na.rm=TRUE), "\n")
  }

  res=list(frame=frame, frame.s=frame.s)
  invisible(res)

}

