
ncs = length(oldprops)
oslopes = matrix(nrow = nrow(SpeciesMat), ncol = ncs)
oslopes = as.data.frame(oslopes)
colnames(oslopes) = c(paste0("OldSlope",oldprops))
oslopes$m0 = oslopes$M = NA

yslopes = matrix(nrow = nrow(SpeciesMat), ncol = length(props))
yslopes = as.data.frame(yslopes)
colnames(yslopes) = c(paste0("YoungSlope",props))
yslopes$m0 = yslopes$M = NA

tslopes = matrix(nrow = nrow(SpeciesMat), ncol = 5)

######## correlations and relative age SD
allprops = paste0(c(rep("_y",length(props)), rep("_o",length(oldprops))),c(props,oldprops))

nc1 = 4*(length(props)+length(oldprops))
corrmat = matrix(nrow = nrow(SpeciesMat), ncol = nc1)
corrmat = as.data.frame(corrmat)
colnames(corrmat) = c(paste0("cor",allprops),paste0("sd_rage",allprops), 
                      as.vector(outer(c("rageL","rageU"),allprops, paste0)))


comparecors = NULL
pdf(paste0(outfolder,"/states/Transform_",name1,len1,"_Strata_",vnum,".pdf"), onefile = T)
par(mfrow=c(2,2))
k = 1
while (k <= nrow(SpeciesMat) ) {
  if(SpeciesMat$RemoveFei[k] == 2) {
    
    k=k+1
    next
  }
  
  spec = SpeciesMat$SpeciesLatinName[k]
  tissue = SpeciesMat$Tissue[k]
  #print(spec)
  if(tissue == "Blood&Skin") idx1 = which(dat1$SpeciesLatinName == spec)
  else if (substr(spec,1,2)=="1.") idx1 = which(dat1$SubOrder == spec & dat1$Tissue == tissue)
  else idx1 = which(dat1$SpeciesLatinName == spec & dat1$Tissue == tissue)
  age1 = dat1$Age[idx1]
  maxage = unique(dat1$maxAgeCaesar[idx1])[1]
  maturity = unique(dat1$AvgMaturity[idx1])[1]
  gt =  unique(dat1$Gestation[idx1])[1]
  
  if(length(unique(age1))==1) {
    SpeciesMat$RemoveFei[k] = 1
    k=k+1
    next
  }
  
  {
    p1 = pmatch(dat1$Basename[idx1], colnames(dat0sesame))
    
    cgidx = pmatch(cgid, rownames(dat0sesame))
    cg1 = dat0sesame[cgidx,p1]
    
  }
  cgmean = colMeans(cg1)
  sdcg = sd(cgmean)
  yslopes$m0[k] = oslopes$m0[k] = m0 = min(cgmean)
  yslopes$M[k] = oslopes$M[k] = M = max(cgmean)   
  
  if(permute) age1 = sample(age1, length(age1), replace = FALSE)
  
  
  slopes = calslope2(cg1,age1, UseMaturity = UseMaturity,
                     maxage= maxage, maturity = maturity,
                    props = props,cut1 = cut1)
  yslopes[k,1:length(props)] = slopes$slopes
  
  res = caloldslope2(cg1,age1, UseMaturity = UseMaturity,
                     maxage = maxage, maturity = maturity,
                    props = oldprops,cut1 = cut1)
  oslopes[k,1:(ncs)] = res$slopes
  
  corrmat[k,] = c(slopes$cors, res$cors,
                  slopes$ragesds, res$ragesds,
                  as.vector(t(slopes$rageRange)),
                  as.vector(t(res$rageRange)))
  
  tage = logli(age1+gt, m1=0.1*maxage+gt)
  tslope1 = calslope2(cg1, tage, UseMaturity = UseMaturity,
                     maxage= maxage, maturity = maturity,
                     props = 0,cut1 = cut1)
  tslopes[k,2] = tslope1$slopes
  tslope2 = caloldslope2(cg1, tage, UseMaturity = UseMaturity,
                     maxage = maxage, maturity = maturity,
                     props = 0,cut1 = cut1)
  tslopes[k,3] = tslope2$slopes
  
  if(length(idx1)>= freq){
    xlab = "Age"
    verboseScatterplot(age1, cgmean,cex.main=cex.main,
                       main=paste0(k,". ",spec,"_", tissue,"\n"),
                       xlab=xlab, ylab= "Mean methylation", type="n")
    text(age1, cgmean, labels = dat1$MammalNumberHorvath[idx1], col = dat1$col.tissue[idx1])
    abline(coef(lm(cgmean~age1)),lty=2,lwd=2)
    
    tslopes[k,4] = cor(age1, cgmean, use = "pairwise.complete")
    # tage = -log(-log((age1+gt)/(maxage+gt)/1.01))
    tage = logli(age1+gt, m1=0.1*maxage+gt)
    
    cor1 = cor(tage, cgmean, use = "pairwise.complete")
    coef1 = signif(coef(lm(cgmean ~ tage)),2)
    plot(tage, cgmean,main = paste0(coef1[1]," / ",coef1[2],"\nCor=",signif(cor1,2)))
    
    coef1 = coef(lm(scale(cgmean) ~ tage))
    tslopes[k,1] = coef1[2]
    tslopes[k,5] = cor1
  }
  k = k+1
}
dev.off()

matslopes = cbind(tslopes, yslopes[,1:length(props)], oslopes[,1:ncs])
colnames(matslopes)[1] = "TSlope"
dim(matslopes)
colSums(is.na(matslopes))
colnames(matslopes) = paste0(name1,len1,colnames(matslopes))
SpeciesMat = cbind(SpeciesMat, matslopes)

colnames(tslopes) = c("TSlope", "TSlope_Young", "TSlope_Old", "Cor_Age", "Cor_TAge")
matslopes = cbind(SpeciesMat[,1:16],tslopes,yslopes,oslopes,corrmat)
write.csv(matslopes, paste0(outfolder,"/states/TSlopes_",name1,len1,vnum,".csv"))


