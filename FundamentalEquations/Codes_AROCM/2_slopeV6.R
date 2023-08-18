library(colorspace)
require(RNOmni)
require(glmnet)
require(caret)
library(plyr)
library(RColorBrewer)
# library(WGCNA)
library(tidyverse)
library(readxl)

outfolder = 'out_AROCM_0.1L'
freq = 5 ##### to include lemurs
anAge <- read.csv("/Users/feiz/Dropbox/MammalianMethCombined/StuffCaesar/anAgeUpdatedCaesarVersion51.csv")
source("Codes_AROCM/1_prepare_data.R")
colSums(is.na(dat1))
dim(dat1)
dim(SpeciesMat)
colSums(is.na(SpeciesMat))
source("lemur.R")

sort(unique(substr(dat1$Folder,1,4)))
tmpd = dat1%>%filter(substr(dat1$Folder,1,4) %in% c("N109"))
table(tmpd$SpeciesLatinName)
# View(dat1%>%filter(substr(dat1$Folder,1,4) %in% c("N119")))
# View(AllSamp%>%filter(substr(AllSamp$Folder,1,4) == "N110"))
## N93, Tursiops aduncus
match("205600840063_R01C01", colnames(dat0sesame))

shrew = dat1%>%filter(SpeciesLatinName=="Sorex cinereus")
table(shrew$Tissue)

shrewnew = rbind(shrew%>%filter(Tissue == "Fetus"),
                 shrew%>%filter(Tissue == "Fetus"))
shrewnew$Tissue = rep(c("Liver","Tail"), each=3)

dat1 = rbind(dat1,shrewnew)
tmp1 = unique(substr(dat1$MammalNumberHorvath,1,3))
tmp1 = strsplit(dat1$MammalNumberHorvath,".",fixed = TRUE)
tmp2 = sapply(tmp1, function(l){
  paste(l[1],l[2], sep=".")
})
sort(unique(tmp2))
dat1$SubOrder = tmp2

write.csv(dat1, paste0(outfolder, "/datSample_slope_",Sys.Date(),".csv"))
########## calculate slopes


permute = FALSE
scaleCpG = TRUE
name1 = ifelse(scaleCpG, "IdentityScaled", "Identity")

states
statesPRC2

library(dplyr)
props = c(1:10)/10  ### 1.00 0.50 0.40 0.30 0.20 0.15 0.10
oldprops = c(1, 1.5, 2)
dim(SpeciesMat)
head(colnames(SpeciesMat), 22)
SpeciesMat = SpeciesMat[,c(1:16)]

source("Codes_AROCM/0_fns_v2.R")
vnum = "v13"
cex.main=1.5
cut1 = freq = 3
j=2
for (j in 1:length(statesPRC2)) {
  cgid = cg_list[[j]]
  len1 = statesPRC2[j]
  
  ## source("Codes_AROCM/3_fitslopeV5.R")
  source("Codes_AROCM/3_fitslopeV8.R")
}

dim(SpeciesMat)
colSums(!is.na(SpeciesMat))[1:42]
####RemoveFei
{
  SpeciesMat$RemoveFei[SpeciesMat$Freq < cut1] = 1
  
  SpeciesMat$RemoveFei[SpeciesMat$MammalNumberHorvath %in% 
                         c("1.3.3","1.3.9", "1.7.1","4.13.2","4.13.11")] = 1
  
  SpeciesMat$RemoveFei[SpeciesMat$MammalNumberHorvath %in% c("1.4.3","6.1.1") &
                         SpeciesMat$Tissue %in% c("Muscle")] = 1
  
  SpeciesMat$RemoveFei[SpeciesMat$MammalNumberHorvath=="4.19.1" &
                         SpeciesMat$Tissue %in% c("Skin")] = 1
  
  SpeciesMat$RemoveFei[SpeciesMat$MammalNumberHorvath=="9.9.1" &
                         SpeciesMat$Tissue %in% c("Brain")] = 1
  SpeciesMat$RemoveFei[SpeciesMat$MammalNumberHorvath=="9.9.3" &
                         SpeciesMat$Tissue %in% 
                         c("Pituitary","Hippocampus", "Hypothalamus")] = 1
  
}

write.csv(SpeciesMat, 
          paste0(outfolder,"/ALL_states_",name1,"_SpeciesTissue_13Slopes_",vnum,".csv"))

### Supplement Table 1
SuppTab1 = read.csv("out0314_TissueSlopes/ALL_states_IdentityScaled_SpeciesTissue_13Slopes_v10.csv",
                    row.names = 1)
head(colnames(SuppTab1),22)
table(SuppTab1$RemoveFei)
SuppTab1 = SuppTab1[SuppTab1$RemoveFei== 0,]
SuppTab1 = SuppTab1[, -(16:17)]
SuppTab1$SpeciesCommonName = NA
m1 = match(SuppTab1$SpeciesLatinName, dat1$SpeciesLatinName)
table(is.na(m1))
SuppTab1$SpeciesCommonName[!is.na(m1)] = dat1$SpeciesCommonName[m1[!is.na(m1)]]
SuppTab1 = SuppTab1[,c(1,718,2:717)]

m2 = match(SuppTab1$SpeciesLatinName, substr(anAgeUse$MammalNumberHorvath,1,3))
SuppTab1$SpeciesLatinName[!is.na(m2)] = anAgeUse$Family[m2[!is.na(m2)]]
View(SuppTab1%>%filter(is.na(SpeciesCommonName)))

SuppTab1$SpeciesCommonName[is.na(SuppTab1$SpeciesCommonName)] = 
  SuppTab1$SpeciesLatinName[is.na(SuppTab1$SpeciesCommonName)]
idx1 = with(SuppTab1,which(SpeciesCommonName == SpeciesLatinName))
SuppTab1$MammalNumberHorvath[idx1] = 
  substr(SuppTab1$MammalNumberHorvath[idx1],1,3)

View(SuppTab1%>%filter(SpeciesCommonName == SpeciesLatinName))

tmpd = read.csv("/Users/feiz/Dropbox/HorvathLabCoreMembers/Josh/ProjectSlope/Article/SupplementalTables/SupplementTable1_v1.csv",
                row.names = 1,check.names=FALSE)
View(cbind(colnames(SuppTab1)[1:16],
           colnames(tmpd)[1:16]))
colnames(SuppTab1)[1:16] = colnames(tmpd)[1:16]
colnames(SuppTab1)[-(1:16)] = substring(colnames(SuppTab1)[-(1:16)],15)

write.csv(SuppTab1%>%
            arrange(SpeciesLatinName), "SupplementTable1_v2.csv",row.names = FALSE)


source("slopes_mat.R")

colnames(corrmat)
dim(corrmat)
ragemat = cbind(SpeciesMat[,1:16], corrmat[,-(1:13)])
colnames(ragemat)
sd_ragemat = ragemat%>%select(contains("sd_rage"))

write.csv(ragemat, paste0(outfolder,"/SpeciesTissueRelativeAge.csv"))

### NEXT
### 4_plots.R


