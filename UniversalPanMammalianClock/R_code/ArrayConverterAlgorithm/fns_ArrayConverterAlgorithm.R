# library(vioplot)
library(glmnet)
library(WGCNA)
library(foreach)
library(doParallel)
library(tidyverse)
#out1 = ArrayConverter_EPICtoMM(cgname,
#                               InputEPIC = betas_epic)
#
ArrayConverter_EPICtoMM <- function(MammalianList = NULL,
                                    InputEPIC, mval = TRUE,
                                    transf = function(a){a},
                                    invfn = function(a){a}){
  nsamp = ncol(InputEPIC)
  
  if(mval) {
    transf = mvalue
    invfn = inv_m
    
  }
  ImpMat = NULL
  for (cgname in MammalianList){
    tmp1 = match(cgname, names(fitall))
    coefmat = fitall[[tmp1]]
    idx1 = match(as.character(coefmat$EPIC_CG[-1]),
                 rownames(InputEPIC))
    if(sum(is.na(idx1)) >0 ) {
      cat(cgname,"Not all EPIC probes for imputation are available.")
      cg_pred = rep(NA, nsamp)
    }
    else {
      xmat = as.matrix(InputEPIC[idx1,])
      if(ncol(xmat)==1){xmat=t(xmat)}
      #print(cgname)
      #print(dim(xmat))
      cg_pred = coefmat$coef[1] + t(xmat)%*%coefmat$coef[-1]
      if(mval) cg_pred = invfn(cg_pred) 
    }
    ImpMat = cbind(ImpMat, cg_pred)
  }
  colnames(ImpMat) = MammalianList
  rownames(ImpMat) = colnames(InputEPIC)
  
  return(ImpMat)
}

mvalue = Vectorize (function(a){
  log(a/(1-a), 2)
}, vectorize.args = c("a"))

inv_m = Vectorize (function(a){
  2^a/(2^a + 1)
}, vectorize.args = c("a"))


corscreen = function(cor1, pick=10000, thr = 1e-3){
  corsort = sort(cor1)
  lower = corsort[pick]
  upper = rev(corsort)[pick]
  zeros = which(abs(cor1) < thr)
  while (length(zeros)< pick) {
    thr = 2*thr
    zeros = which(abs(cor1) < thr)
  }
  zeros = sample(zeros, pick)
  
  idxpick = sort(unique(c(which(cor1<= lower), which(cor1>= upper), zeros)))
  return(idxpick)
}


ImpMM_v1 <- function(cgname, idx1, train.idx, test.idx,
                     width = 1e6, geneRegion = FALSE,
                     doplot=FALSE){
  
  cg1 = as.numeric(normalized_betas_sesame[rownames(normalized_betas_sesame) == cgname,idx1])
  
  chr = as.character(mm_annot$seqnames[mm_annot$CGid == cgname])
  
  if(geneRegion){
    geneStart = mm_annot$geneStart[mm_annot$CGid == cgname]
    geneEnd = mm_annot$geneEnd[mm_annot$CGid == cgname]
    
    cgs_epic = epic_annot$CGid[epic_annot$seqnames == chr &
                                 epic_annot$start >= geneStart &
                                 epic_annot$start <= geneEnd] 
    
  }else {
    CGstart = mm_annot$CGstart[mm_annot$CGid == cgname]
    cgs_epic = epic_annot$CGid[epic_annot$seqnames == chr &
                                 epic_annot$start >= CGstart - width &
                                 epic_annot$start <= CGstart + width] 
    
  }
  
  cgs_epic = intersect(cgs_epic, rownames(betas_epic))
  npred = length(cgs_epic)
  if(npred<2) return(c(npred, rep(median(cg1[train.idx]),length(idx1))))
  
  epic_mat = betas_epic[cgs_epic, idx1]
  # verboseScatterplot(cg1, betas_epic[cgname, idx1])
  
  m1 = cv.glmnet(t(epic_mat[,train.idx]), cg1[train.idx], alpha=.5)
  # plot(m1)
  # coef(m1,s="lambda.min")
  
  cg_fitted = predict(m1, newx = t(epic_mat), s="lambda.min" )  
  
  if(doplot){
    # verboseScatterplot(cg1, cg_fitted)
    verboseScatterplot(cg1[test.idx], cg_fitted[test.idx])
    
  }
  return(c(npred,cg_fitted))
}


ImpMM_v4 <- function(cgname, idx1, train.idx, test.idx,
                     width = 1e6, geneRegion = FALSE, use30k = FALSE,
                     transf = function(a){a}, invfn = function(a){a},
                     doplot=FALSE){
  
  cg1 = as.numeric(normalized_betas_sesame[rownames(normalized_betas_sesame) == cgname,idx1])
  cg1 = transf(cg1)
    
  chr = as.character(mm_annot$seqnames[mm_annot$CGid == cgname])
  
  if(geneRegion){
    geneStart = mm_annot$geneStart[mm_annot$CGid == cgname]
    geneEnd = mm_annot$geneEnd[mm_annot$CGid == cgname]
    
    cgs_epic = epic_annot$CGid[epic_annot$seqnames == chr &
                                 epic_annot$start >= geneStart &
                                 epic_annot$start <= geneEnd] 
    
  }else if (use30k) {
    cor1 = cor(cg1[train.idx], t(betas_epic[, idx1[train.idx]]))
    idxpick = corscreen(cor1)
    cgs_epic = rownames(betas_epic)[idxpick]
    
  } else {
    CGstart = mm_annot$CGstart[mm_annot$CGid == cgname]
    cgs_epic = epic_annot$CGid[epic_annot$seqnames == chr &
                                 epic_annot$start >= CGstart - width &
                                 epic_annot$start <= CGstart + width] 
    
  }
  
  cgs_epic = intersect(cgs_epic, rownames(betas_epic))
  npred = length(cgs_epic)
  if(npred<2) {
    y = invfn(median(cg1[train.idx]))
    return(c(npred, rep(y,length(idx1)) ))
  }
  
  epic_mat = betas_epic[cgs_epic, idx1]
  # verboseScatterplot(cg1, betas_epic[cgname, idx1])
  
  m1 = cv.glmnet(t(epic_mat[,train.idx]), cg1[train.idx], alpha=.5)
  # plot(m1)
  # coef(m1,s="lambda.min")
  
  cg_fitted = predict(m1, newx = t(epic_mat), s="lambda.min" )  
  
  if(doplot){
    # verboseScatterplot(cg1, cg_fitted)
    verboseScatterplot(cg1[test.idx], cg_fitted[test.idx])
    
  }
  
  y = invfn(cg_fitted)
  return(c(npred,y))
}

fitMM_v1 = function(cgname, idx1, train.idx, test.idx=NULL,
                     width = 6e7, uselambda = "lambda.min",
                     # transf = function(a){a}, invfn = function(a){a},
                     doplot=FALSE, returnCoef = TRUE){
  mvalue = Vectorize (function(a){
    log(a/(1-a), 2)
  }, vectorize.args = c("a"))
  
  inv_m = Vectorize (function(a){
    2^a/(2^a + 1)
  }, vectorize.args = c("a"))
  transf = mvalue
  invfn = inv_m
  
  cg1 = as.numeric(normalized_betas_sesame[rownames(normalized_betas_sesame) == cgname,idx1])
  cg1 = transf(cg1)
  
  chr = as.character(mm_annot$seqnames[mm_annot$CGid == cgname])
  
  {
    CGstart = mm_annot$CGstart[mm_annot$CGid == cgname]
    cgs_epic = epic_annot$CGid[epic_annot$seqnames == chr &
                                 epic_annot$start >= CGstart - width &
                                 epic_annot$start <= CGstart + width] 
    
  }
  
  cgs_epic = intersect(cgs_epic, rownames(betas_epic))
  npred = length(cgs_epic)
  if(npred<2) {
    cat("Not enough EPIC CGs for imputation, try longer bandwidth.")
    return(c(npred ))
  }
  
  epic_mat = t(betas_epic[cgs_epic, idx1])
  m1 = cv.glmnet(epic_mat[train.idx,], cg1[train.idx], alpha=.5)
  if(doplot){
    plot(m1)
  }
  coef1 = coef(m1, s=uselambda)
  coefs = data.frame(EPIC_CG = rownames(coef1)[coef1[,1]!=0],
             coef = coef1[coef1!=0])
  cg_fitted = predict(m1, newx = epic_mat, s=uselambda )  
  y = invfn(cg_fitted)
  
  if(returnCoef) return(coefs)
  else return(y)
}

imputeMM_v1 = function(cglist, 
                       doparallel = FALSE){
  
}
