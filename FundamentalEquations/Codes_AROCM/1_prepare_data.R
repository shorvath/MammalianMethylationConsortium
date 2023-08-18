

# source("../Data/fns.R")
confidence = FALSE
# source("../Data/readdata_run.R")
AllSamp = read.csv("/Users/feiz/Dropbox/MammalianMethCombined/StuffCaesar/datSampleCombinedfinal.csv")

tmpdat = anAge %>% 
  select(SpeciesLatinName,Order,
         maxAgeCaesar,weightCaesar,Gestation.Incubation..days.,
         averagedMaturity.yrs,MammalNumberHorvath,
         (contains("litter") & !contains("PanTHERIA")))
colnames(tmpdat)

AllSamp = merge(AllSamp,tmpdat, by.x = "SpeciesLatinName")
AllSamp$Gestation = AllSamp$Gestation.Incubation..days./365
AllSamp$AvgMaturity = AllSamp$averagedMaturity.yrs
table(AllSamp$SpeciesLatinName[is.na(AllSamp$AvgMaturity)])
table(AllSamp$Folder[is.na(AllSamp$AvgMaturity)])

traningIndex <- read.csv("/Users/feiz/Dropbox/HorvathLabCoreMembers/Ake/UniversalClock/TrainingIndex_final.csv")
colnames(traningIndex)


if(FALSE){
  datN94 = AllSamp%>%filter(substr(Folder,1,3) == "N94")%>%
    select(OriginalOrderInBatch,Basename, Age,ConfidenceInAgeEstimate,
           CanBeUsedForAgingStudies,Tissue)%>%
    mutate(tissueColor = unique(traningIndex$tissueColor[traningIndex$Tissue=="Blood"]),
           SmallsizeSpecies = 0,
           Training = NA,
           LOFO.id = NA)
  
  
  traningIndex <- rbind(traningIndex, datN94)
  
}

############ combined with marsupials
mars = read.csv("/Users/feiz/Dropbox/HorvathLabCoreMembers/Ake/UniversalClock/TrainingIndex_final_MarsupialsOnly.csv")
mars = mars%>%select(colnames(traningIndex))
identical(colnames(mars), colnames(traningIndex))
traningIndex = rbind(traningIndex, mars)
colnames(traningIndex)

table(traningIndex$SmallsizeSpecies)
tmp1 <- c(pmatch(traningIndex$Basename,AllSamp$Basename))
summary(tmp1)
dat1 = AllSamp[tmp1,]
dat1$col.tissue <- traningIndex$tissueColor
dat1$Training <- traningIndex$Training
dat1$LOFO.id <- traningIndex$LOFO.id
dat1$SmallsizeSpecies = traningIndex$SmallsizeSpecies
colSums(is.na(dat1))
sort(table(dat1$Order))

dat1 = dat1%>%filter(ConfidenceInAgeEstimate>= 90)
dat1$Tissue[dat1$Tissue=="Spleen"] = "Blood"
dat1$Tissue[dat1$Tissue=="Ear"] = "Skin"

load("/Users/feiz/Dropbox/MammalianMethCombined/StuffCaesar/NormalizedData/all_probes_all_samples_sesame.RData")
if (!is.numeric(dat0sesame[,1])) {
  rownames(dat0sesame) <- dat0sesame[,1]
  dat0sesame <- dat0sesame[,-1]
} 

id1 = substr(rownames(dat0sesame),1,2) != "rs"
table(id1)
dat0sesame = dat0sesame[id1,]

tail(rownames(dat0sesame))
dim(dat0sesame)
##### species tissue mat

{
  ### input: cgid, dat1

  print(freq)
  # dat1 = AllSamp
  tab1 = data.frame(with(dat1,table(SpeciesLatinName,Tissue)))
  tab1 = tab1[tab1$Freq>= freq,]
  #colnames(tab1)[1] = "SpeciesLatinName"
  SpeciesMat <- ddply(dat1,c("SpeciesLatinName","Tissue"),function(dat1){
    if(nrow(dat1)>= freq & length(unique(dat1$Age))>1   ){
      tmp1 = dat1%>% select(MammalNumberHorvath,
                            col.tissue,Order,
                            contains("litter"),
                            weightCaesar,Gestation,
                            AvgMaturity,maxAgeCaesar)
      tmp2 = range(dat1$Age)/unique(tmp1$maxAgeCaesar)
      data.frame(tmp1[1,],ageRangeL = tmp2[1],ageRangeU = tmp2[2])
    }
  })
  dim(SpeciesMat)
  # View(SpeciesMat%>%filter(SpeciesLatinName=="Sorex cinereus"))
  SpeciesMat<- merge(SpeciesMat, tab1, by=c("SpeciesLatinName","Tissue"))
  SpeciesMat$MammalNumberHorvath = as.character(SpeciesMat$MammalNumberHorvath)
  SpeciesMat$SpeciesLatinName = as.character(SpeciesMat$SpeciesLatinName)
  SpeciesMat$Tissue = as.character(SpeciesMat$Tissue)
  
  
  # SpeciesMat[is.na(SpeciesMat$Litters.Clutches.per.year),]
  # SpeciesMat[is.na(SpeciesMat$Litter.Clutch.size),]
  SpeciesMat$Litter.sizeXclutchperyear = SpeciesMat$Litter.Clutch.size*SpeciesMat$Litters.Clutches.per.year
  
  
}
colnames(SpeciesMat)
summary(SpeciesMat%>%select(contains("litter")))
# SpeciesMat%>%filter(is.na(Litters.Clutches.per.year))

table(SpeciesMat$ageRangeL == SpeciesMat$ageRangeU)
head(sort(SpeciesMat$ageRangeU),11)
SpeciesMat$RemoveFei = 0
SpeciesMat$RemoveFei[SpeciesMat$ageRangeU< 0.1] = 1
# View(SpeciesMat%>%filter(RemoveFei == 1))
# SpeciesMat%>% filter(SpeciesLatinName == 'Mus musculus' & Tissue == 'Fibroblast')
#SpeciesMat$RemoveFei[SpeciesMat$SpeciesLatinName == 'Mus musculus' & 
#                      SpeciesMat$Tissue == 'Fibroblast'] = 1

SpeciesMat <- SpeciesMat[order(SpeciesMat$maxAgeCaesar),]
write.csv(SpeciesMat, paste0(outfolder,"/SpeciesTissueStrata.csv"))

######### marsupial arrays
if(FALSE){
  N26dat0custom = readRDS("/Users/feiz/Dropbox/N26.2019-9088TasmanianDevilBloodCarolynHogg/NormalizedDataCustomSesame/all_probes_sesame_normalized.RDS")
  N42dat0custom = readRDS("/Users/feiz/Dropbox/N42.2019-9211OpossumKarenSears/NormalizedDataCustomSesame/all_probes_sesame_normalized.RDS")
  load("/Users/feiz/Dropbox/N81.MacropusFromN27/NormalizedDataCustomSesame/all_probes_sesame_normalized.Rdata")
  N81dat0custom = normalized_betas_sesame
  p1 = pmatch(colnames(N81dat0custom)[-1], dat1$Basename)
  p1 = p1[!is.na(p1)]
  table(dat1$SpeciesLatinName[p1])
  
  identical(rownames(N26dat0custom), rownames(N42dat0custom) )
  identical(rownames(N26dat0custom), rownames(N81dat0custom) )
  dat0custom = cbind(N26dat0custom, N42dat0custom[,-1])
  
  saveRDS(dat0custom,"dat0custom_marsupials.rds")
  
  cglist_amin = read.csv("/Users/feiz/Dropbox/HorvathLabCoreMembers/Josh/ProjectSlope/stackHMM status of Eutherian CpGs.AminV1.csv",row.names = 1)
  cglist_amin2 = read.csv("stackHMM status of Eutherian CpGs.Amin.csv",row.names = 1)
  identical(cglist_amin$CGid, cglist_amin2$CGid)
  
  PRCbound = read.csv("/Users/feiz/Dropbox/MammalianArrayNormalizationTools/ChromatinStatesJasonErnstSooBinKwon/StackHMM, repeat elements, PRCchipseq/mammalian_hg19_fullStack_repeats_PRCchipseq, aggregated.csv")
  head(PRCbound)
  tmptab = with(PRCbound,table(full_stacked_state,PRC1_Amin))
  tmptab[tmptab[,1]>100,]
  tmptab = with(PRCbound,table(full_stacked_state,PRC2_Amin))
  tmptab[tmptab[,1]>100,]
  head(with(PRCbound,table(full_stacked_state,PRC2_Amin)), 22)
  
  head(cglist_amin)
  table(cglist_amin$CpGislandIn48PlusSpecies[cglist_amin$stackHMM=="ReprPC1"])
  p1 = match(cglist_amin$CGid, PRCbound$CGid)
  summary(p1)
  CGstates = merge(cglist_amin,PRCbound[p1,-(1:4)],by = 'CGid')
  table(CGstates$stackHMM == CGstates$full_stacked_state)
  View(CGstates%>% filter(stackHMM !=full_stacked_state))
  table(CGstates$PRC2_Amin[CGstates$full_stacked_state=="ReprPC1"])
  
}




############## cg lists by states
# BivProm1, BivProm2, PromF5, TSS2, TxEx4
# states = rev(names(tab1)[tab1>= numcgs])
CGstates = read.csv("/Users/feiz/Dropbox/HorvathLabCoreMembers/Josh/ProjectSlope/StuffFei/mammalian_hg19_fullStack_repeats_PRCchipseq_updated.csv",
                    row.names = 1)
head(CGstates)

tab1 = sort(table(CGstates$full_stacked_state))
numcgs = 5
states = statesPRC2 = NULL
cg_list = list()
jj=0
for (j in length(tab1):1 ) {
  cs = names(tab1)[j]
  
  cgid = CGstates$CGid[CGstates$full_stacked_state==cs & CGstates$PRC2_Amin == 1]
  if(length(cgid)>= numcgs){
    jj=jj+1
    cg_list[[jj]] = cgid
    states = c(states,cs)
    statesPRC2 = c(statesPRC2,paste0(cs,"+"))
  }
  cgid2 = CGstates$CGid[CGstates$full_stacked_state==cs & CGstates$PRC2_Amin == 0]
  if(length(cgid2)>= 100){
    jj=jj+1
    cg_list[[jj]] = cgid2
    states = c(states,cs)
    statesPRC2 = c(statesPRC2,paste0(cs,"-"))
  }
  
}
cbind(states, statesPRC2)
length(cg_list )
sapply(cg_list, length)
# setdiff(states, NegStates)

statesTab = data.frame(states = states, statesPRC2=statesPRC2,
                       Freq = sapply(cg_list, length))

write.csv(statesTab, paste0(outfolder,"/ChromStates.csv"))
saveRDS(cg_list, paste0(outfolder,"/ChromStates_CGlist.rds"))
