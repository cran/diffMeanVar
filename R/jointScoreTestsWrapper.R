# revised on Dec. 12, 2015
#  (1) replace 'y' by 'group'
#  (2) replace 'x' by 'group'
#  (3) rename 'jointScoreTests' to 'jointScoreTestsWrapper'
#  (4) rename 'jointScoreTests.default' to 'jointScoreTests'
#
# created on July 17, 2015
#  (1) joint test for mean or variance difference by
#      using score test of logistic regression 
#         to test mean or variance difference 

# group - vector of binary values
# value - continuous variable
jointScoreTests= function(
  value, 
  group, 
  meanTestFunc = scoreTestMean,
  varTestFunc = AWvarTest)
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

  ########################
  # test statistic for mean difference
  ########################
  res.mean=meanTestFunc(value, group)
  U1=res.mean$U1
  varU1=res.mean$varU1
  x2=res.mean$x2
  x2bar=res.mean$x2bar

  ########################
  # test statistic for variance difference
  ########################
  res.var=varTestFunc(value, group)
  U2=res.var$U2
  varU2=res.var$varU2
  z=res.var$z
  zbar=res.var$zbar

  ###################
  ##############
  groupbar = mean(group, na.rm=TRUE)
  covU1U2 = groupbar*(1-groupbar)*sum((x2-x2bar)*(z-zbar), na.rm=TRUE)

  #############
  # score test statistic
  #############
  sigma11 = varU1
  sigma12 = covU1U2
  sigma22 = varU2

  denom = sigma11*sigma22 - sigma12^2

  sigma11.inv = sigma22/denom
  sigma12.inv = -sigma12/denom
  sigma22.inv = sigma11/denom

  Tjoint = U1^2*sigma11.inv + 2*U1*U2*sigma12.inv + U2^2*sigma22.inv

  pval=1-pchisq(Tjoint, df=2)

  res=list(U1=U1, varU1=varU1, U2=U2, varU2=varU2, covU1U2=covU1U2,
    stat=Tjoint, pval=pval)
  return(res)

}

# esFlag = "es" (for ExpressionSet object) or "meth" (for methylumi object)
jointScoreTestsWrapper=function(
  es,
  grpVar = "group", 
  meanTestFunc = scoreTestMean,
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
    ans=try(resi <- jointScoreTests(value=xi,group=y, 
      meanTestFunc=meanTestFunc,
      varTestFunc=varTestFunc))
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

