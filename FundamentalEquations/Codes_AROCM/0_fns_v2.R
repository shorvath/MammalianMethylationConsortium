logplus = function(y){
  min1 = min(y[y>0 & !is.na(y)])/2
  y[y<=0 & !is.na(y)] = min1
  return(y)
}

verboseScatterplot = function(x,y,main="",
         xlim=range(x,na.rm = TRUE),ylim=range(y,na.rm = TRUE),
         ...){
  cor1 = cor.test(x,y, use = "pairwise.complete")
  
  p1 = cor1$p.value
  cor2 = round(cor1$estimate,2)
  
  plot(x,y, xlab=xlab,ylab=ylab,pch=19,
       xlim=xlim,ylim=ylim,
       main=paste0(main,"Cor=",cor2,
                        "\nP=",signif(p1,3)))
}

invfn <- function(x){
  m = mean(x, na.rm=TRUE)
  mae = mean(abs(x-m), na.rm=TRUE)
  medae = median(abs(x-m), na.rm=TRUE)
  sdx = sd(x, na.rm=TRUE)
  return(c(mae,medae,sdx))
}

caloldslope2 <- function(cg1,age1, UseMaturity = TRUE,
                         maxage=100,maturity, props = c(1, 1.5),
                        cut1=5,intercept = TRUE, calslope=TRUE){
  if(ncol(cg1) != length(age1)) return(print("Input lengths differ."))
  n = ncol(cg1)
  # y = cg1
  
  slopes = cors = SD = SDcg = rep(NA, length(props))
  rageRange = matrix(nrow = length(props), ncol = 2) 
  for(j in 1:length(props)){
    prop = props[j]
    if (UseMaturity) id1 = which(age1>= prop*maturity) else
      id1 = which(age1>= prop*0.1*maxage)
    n1 = length(id1)
    
    if (n1>= cut1) {
      sdcg = sd(colMeans(cg1[,id1]))
      y = scale(colMeans(cg1[,id1]))
      slope1 = ifelse(intercept, 
                      coef(lm(y~age1[id1]))[2], 
                      coef(lm(y~age1[id1]-1)) ) 
      slopes[j] = slope1
      SD[j] = sd1 = sd(age1[id1]/maxage)
      SDcg[j] = sdcg
      cor1 = cor(y, age1[id1])
      cors[j] = cor1
      
      rageRange[j,] = range(age1[id1])/maxage
    }
  }
  return(list(slopes = slopes, cors = cors,
              ragesds = SD,
              SDMeth = SDcg,
              rageRange = rageRange))
}

calslope2 <- function(cg1,age1, UseMaturity = FALSE,
                      maxage=100,maturity, props = props,
                      cut1=5, intercept = TRUE){
  if(ncol(cg1) != length(age1)) return(print("Input lengths differ."))
  n = ncol(cg1)
  # y = cg1
  
  slopes = cors = SD = SDcg = rep(NA, length(props))
  rageRange = matrix(nrow = length(props), ncol = 2) 
  for(j in 1:length(props)){
    prop = props[j]
    if(UseMaturity& j<= 3){
      id1 = which(age1<= 10*prop*maturity)
      n1 = length(id1)
      
    }else {
      id1 = which(age1<= prop*maxage)
      n1 = length(id1)
      
    }
    
    if (n1>= cut1) {
      sdcg = sd(colMeans(cg1[,id1]))
      y = scale(colMeans(cg1[,id1]))
      
      slope1 = ifelse(intercept, 
                      coef(lm(y~age1[id1]))[2], 
                      coef(lm(y~age1[id1]-1)) ) 
      slopes[j] = slope1
      SD[j] = sd1 = sd(age1[id1]/maxage)
      SDcg[j] = sdcg
      cor1 = cor(y, age1[id1])
      cors[j] = cor1
      
      rageRange[j,] = range(age1[id1])/maxage
    }
  }
  return(list(slopes = slopes, cors = cors,
              ragesds = SD,
              SDMeth = SDcg,
              rageRange = rageRange))
}


AROCM = function(cgmean, age1, intercept = TRUE){
  {
    sdcg = sd(cgmean)
    y = scale(cgmean)
    
    if(intercept) lm1 = lm(y~age1) else
      lm1 = lm(y~age1-1)
    slope1 = ifelse(intercept, 
                    coef(lm1)[2], 
                    coef(lm1) )
    r2 = summary(lm1)$r.squared
    
    cor1 = cor(y, age1)
    
    out = c(slope1, cor1, sd(age1), sdcg, r2)
    names(out) = c("AROCM", "Cor", "SD_Age", "SD_Methyl", "R2")
  }
  return(out)

}
fitAROCM = function(dat1, SpeciesMat, ageRange = c(0,1), 
                    relativeAge = TRUE, plotout = FALSE,
                    cut1 = 3){
  
  # outmat = SpeciesMat[,c(1:3, 9:12)]
  outmat = matrix(NA, nrow = nrow(SpeciesMat), ncol = 5)
  k = 1
  while (k <= nrow(SpeciesMat) ) {
    
    spec = SpeciesMat$SpeciesLatinName[k]
    tissue = SpeciesMat$Tissue[k]
    #print(spec)
    if(tissue == "Blood&Skin") idx1 = which(dat1$SpeciesLatinName == spec)
    else if (substr(spec,1,2)=="1.") idx1 = which(dat1$SubOrder == spec & dat1$Tissue == tissue)
    else idx1 = which(dat1$SpeciesLatinName == spec & dat1$Tissue == tissue)
    age1 = dat1$Age[idx1]
    maxage = unique(dat1$maxAgeCaesar[idx1])[1]
    maturity = unique(dat1$AvgMaturity[idx1])[1]
    
    
    {
      p1 = pmatch(dat1$Basename[idx1], colnames(dat0sesame))
      
      cgidx = pmatch(cgid, rownames(dat0sesame))
      cg1 = dat0sesame[cgidx,p1]
      
    }
    cgmean = colMeans(cg1)
    
    if(relativeAge) 
      idx2 = which(age1>= ageRange[1]*maturity & age1<= ageRange[2]*maxage) else 
        idx2 = which(age1>= ageRange[1] & age1<= ageRange[2])
    
    if(length(unique(age1[idx2]))>= cut1){
      out1 = AROCM(cgmean[idx2], age1[idx2])
      out1[5] = out1[3]/maxage
      names(out1)[5] = "SD_RAge"
      outmat[k,] = out1
    }
    
    k = k+1
  }
  colnames(outmat) = names(out1)
  
}

plotspecs = function(mat1, nc, prop = 1,letter=4,title1=NA,tit1 = NA,
                     zoom = NA, a=NA, b=NA){
  mat1 = mat1[!is.na(mat1[,nc]),]
  x = mat1$maxAgeCaesar
  y = mat1[,nc]
  xvar = "a/Lifespan^b"
  if(prop<1) x = prop*x
  # if(prop>1) x = as.numeric(mat1$AvgMaturity)
  xvar1 = "Max Lifespan"
  # xlab1 = ifelse(prop<1, paste0(prop,xvar1), ifelse(prop==1, paste0(xvar1), 'AvgMaturity' ))
  xlab1 = ifelse(prop<1, paste0(prop,xvar1), paste0(xvar1))
  if(is.na(tit1)){
    tit1 = ifelse(prop<=1, paste0("(L,U)=(0,",prop,")"), 
                  paste0("(L,U)=(1.5*Maturity,Lifespan)"))
  }
  
  tmp1 = strsplit(colnames(mat1)[nc],".",fixed = T)[[1]]
  ylab1 = ifelse(tmp1[2]=="Slope", paste0(tmp1[1], tmp1[2]), paste0(tmp1[2], tmp1[3]) )
  
  
  offset = ifelse(min(y)<0, -min(y)*1.5, 0)
  y = y + offset 
  
  if(is.na(b)){
    lm1 = lm(log(y)~log(x))   #### log(x) = a*log(y) + b
    coef1 = as.numeric(coef(lm1))
    x1 = exp(coef1[1])*(x)^coef1[2]
    a = round(exp(coef1[1]),2)
    b = round(-coef1[2],2)
  }else if (is.na(a)){
    tmpx =  (x)^(-b)
    lm1 = lm(y~ tmpx-1, weights = tmpx)    
    coef1 = as.numeric(coef(lm1))
    x1 = coef1[1]*(x)^(-b)
    a = round(coef1[1],2)
  }else {
    x1 = a*(x)^(-b)
    
  }
  
  # name1 = paste(unlist(strsplit(colnames(mat1)[nc],".", fixed = TRUE))[c(1,3)],collapse = "")
  
  if(is.na(zoom)){
    cor1 = cor(x1,y,use = "pairwise.complete", method = "p")
    cor2 = cor(x1,y,use = "pairwise.complete", method = "s")
    plot(x1,y, type="n",xlab = paste0("a/(",xlab1,")^b"),ylab = paste0('Slope_',title1),
         main = paste0('N=',length(y),", ",tit1,'\na=',a, ', b=',b, 
                       '\n',"sCor=",round(cor2,3),", pCor=",round(cor1,3)) )
    text(x1,y, labels = mat1$MammalNumberHorvath)
    abline(0,1,lty=2)
    abline(coef(lm(y~x1)))
  } else {
    idx1 = which(x1 <= zoom)
    x1 = x1[idx1]
    y = y[idx1]
    cor1 = cor(x1,y,use = "pairwise.complete", method = "p")
    cor2 = cor(x1,y,use = "pairwise.complete", method = "s")
    xr = range(c(x1,y), na.rm = T)
    plot(x1,y, type="n",xlab = xlab1,ylab = 'Slope',xlim=xr,ylim=xr,
         main = paste0(letters[letter],'. ',xvar,"<",zoom,"\nsCor=",round(cor2,3),", pCor=",round(cor1,3)) )
    text(x1,y, labels = mat1$MammalNumberHorvath[idx1])
    abline(0,1,lty=2)
    abline(coef(lm(y~x1)))
  }
  
  return(coef1[1])
}


combineTissues <- function(slopes_mat, idx1, nc){
  speciesSlopes = ddply(slopes_mat[idx1,], c("SpeciesLatinName"), 
                        function(mat1){
                          tmpd = apply(mat1[,nc],2,median, na.rm=T)
                          # Freq = sum(mat1$Freq)
                          # c(Freq, tmpd)
                          tmpd
                        })
  
  tmp1 = match(speciesSlopes$SpeciesLatinName, SpeciesMat$SpeciesLatinName)
  speciesSlopes$MammalNumberHorvath = SpeciesMat$MammalNumberHorvath[tmp1]
  speciesSlopes
}

caloldslope <- function(cg1,age1, maxage=100,maturity, props = c(1, 1.5),
                     intercept = TRUE, calslope=TRUE){
  if(length(cg1) != length(age1)) return(print("Input lengths differ."))
  
  n = length(cg1)
  y = cg1
  
  slopes = cors = SD = rep(NA, length(props))
  rageRange = matrix(nrow = length(props), ncol = 2) 
  for(j in 1:length(props)){
    prop = props[j]
    id1 = which(age1>= prop*maturity)
    n1 = length(id1)
    
    if (n1>= 5) {
      slope1 = ifelse(intercept, 
                      coef(lm(y[id1]~age1[id1]))[2], 
                      coef(lm(y[id1]~age1[id1]-1)) ) 
      slopes[j] = slope1
      SD[j] = sd1 = sd(age1[id1]/maxage)
      cor1 = cor(y[id1], age1[id1])
      cors[j] = cor1
      
      rageRange[j,] = range(age1[id1])/maxage
    }
  }
  return(list(slopes = slopes, cors = cors,
              ragesds = SD,
              rageRange = rageRange))
}

calslope <- function(cg1,age1, maxage=100,maturity, props = props,
                     intercept = TRUE){
  if(length(cg1) != length(age1)) return(print("Input lengths differ."))
  
  n = length(cg1)
  y = cg1
  
  slopes = cors = SD = rep(NA, length(props))
  rageRange = matrix(nrow = length(props), ncol = 2) 
  for(j in 1:length(props)){
    prop = props[j]
    if(prop>1){
      id1 = which(age1>= prop*maturity)
      n1 = length(id1)
      
    }else {
      id1 = which(age1<= prop*maxage)
      n1 = length(id1)
      
    }
    if (n1>= 5) {
      slope1 = ifelse(intercept, 
                      coef(lm(y[id1]~age1[id1]))[2], 
                      coef(lm(y[id1]~age1[id1]-1)) ) 
      slopes[j] = slope1
      SD[j] = sd1 = sd(age1[id1]/maxage)
      cor1 = cor(y[id1], age1[id1])
      cors[j] = cor1
      
      rageRange[j,] = range(age1[id1])/maxage
    }
  }
  return(list(slopes = slopes, cors = cors,
              ragesds = SD,
              rageRange = rageRange))
}
  
slopefn <- function(cg1,age1, FUN = revGrowth,intercept = TRUE, takeLog = FALSE, Hampel=TRUE,
                    c0=0,c1=1, c2 = 1.001,s=4,plus=FALSE, plot = FALSE,plotname="Identity",
                    ylab=len1,
                    cex.main = 1, returnCor = FALSE){
  if(length(unique(age1))==1) return(NA)
  m0 = min(cg1)*c1
  M = ifelse(plus, min(max(cg1)+c2-1,1), min(max(cg1)*c2,1)) 
  if(takeLog) {
    if(min(age1)>0) c0=0
    else if(min(age1)==0) c0=min(0.01, min(age1[age1>0]))
    else c0 = -1.1*min(age1)
    age1 = log(age1 + c0)
    if(c0>0) print(c0)
  }
  
  y = FUN(cg1, m0, M)
  if(Hampel) y = hampel(y,s=s)
  slope1 = ifelse(intercept, coef(lm(y~age1))[2], coef(lm(y~age1-1)) ) 
  coef1 = coef(lm(y~age1))
  cor1 = cor(y,age1)
  
  if(plot){
    xlab = ifelse(takeLog, "LogAge", "Age")
    verboseScatterplot(age1, y,cex.main=cex.main,
                       main=paste0(k,". ",spec,"_", tissue,"\n",plotname,", slope=",signif(slope1,2),"\n"),
                       xlab=xlab, ylab= ylab, type="n")
    text(age1, y, labels = dat1$MammalNumberHorvath[idx1], col = dat1$col.tissue[idx1])
    abline(coef1,lty=2,lwd=2)
  }
  return(ifelse(returnCor,cor1,slope1) )
}

  