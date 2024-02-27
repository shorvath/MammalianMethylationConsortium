
## Since Mammalian lifespan predictor is trained on species-level data, we need to average all-sample data by species
## This is just a singled out convenience step, produce the averaged data and save it,
## so we don't have to run it many times for subsequent different predictor training
## In the end, save the file to an RDS file to be used for other Elastic Net training code.


library(dplyr)


combinedsamples = read.csv("./Samples_used_for_LifespanProject_v2corrected.csv", stringsAsFactors = F)
source('./utilities_analysis.R')


# remember to exclude N65.2020-9118MouseCalico for lifespan project
# and 86 87 88 89 62 for lifespan

useNewCode = T
noExperiments = F

print(paste0("Number of samples after filtering out N70 and hybrids: ", nrow(combinedsamples)))

### filtering no NAs in maxAgeCaesar

anageCaesar = read.csv("~/StuffCaesar/anAgeUpdatedCaesarVersion49.csv", stringsAsFactors = F)
combinedsamples = combinedsamples %>% left_join(anageCaesar, by = "SpeciesLatinName")

## exclude folders that we cannot use (copyright), non-mammals, and filter out iPSC, ES cells, as well as hybrid species like Tt x Tt gilli
combinedsamples = combinedsamples[!combinedsamples$SpeciesLatinName == "", ] 
combinedsamples = combinedsamples[!is.na(combinedsamples$SpeciesLatinName) & !combinedsamples$SpeciesLatinName == "Peromyscus hybrid polionatus maniculatus", ] 

# exclude experiment folders 
## Note: these are the folders that I excluded IN ADDITION TO STEVE'S REMOVE FOLDER
combinedsamples = combinedsamples[!grepl("N65\\.|N70\\.|N74\\.|N84\\.|N86\\.|N87\\.|N88\\.|N89\\.|N62\\.", combinedsamples$Folder), ]
combinedsamples = combinedsamples[!grepl("N02\\.|N25\\.|N39\\.|N112", combinedsamples$Folder), ]
####
# exclude new data and non-human data (TO BE MODIFIED IN THE FUTURE)
combinedsamples = combinedsamples[!grepl("N95|N96|N97|N98|N99|N100|N101|N102|N103|N104|N105|N107|N108|N111|N113", 
                                         combinedsamples$Folder), ]
#
## Add Steve's exclude Dec2021
combinedsamples = combinedsamples[!grepl("N71\\.|N79\\.|N80\\.|N82\\.|N83\\.|N85\\.|N128\\.", combinedsamples$Folder), ]

## ADD only N141
combinedsamples = combinedsamples[!grepl("N138\\.|N139\\.|N140\\.|N142\\.", combinedsamples$Folder), ]

combinedsamples = combinedsamples[!grepl("iPSC|ES", combinedsamples$Tissue), ]
combinedsamples = combinedsamples[combinedsamples$CanBeUsedForAgingStudies == "yes", ]
combinedsamples = combinedsamples[!grepl("Tt x Tt gilli", combinedsamples$SpeciesLatinName), ]

## SwordFish N122, frogs, Frog folders N95 N140
combinedsamples = combinedsamples[!grepl("N95\\.|N122\\.|N140\\.", combinedsamples$Folder), ]


######
## This is just a convenience function to help us avarege normalized CpG data by species
NoobAverageSpecies <- function(samples, dt = NA, subsetVar = c(NA, NA), saveProof = NA) {
  library(dplyr)
  
  if(!is.na(subsetVar[1])) {
    samples = samples[!is.na(samples[, subsetVar[1]]) & samples[, subsetVar[1]] %in% subsetVar[-1], ]
    dat01 <- dat0sesame[rownames(dat0sesame) %in% samples$Basename, ]
  } else {
    dat01 <- dat0sesame 
  }
  print(dim(dat01))
  print(dim(samples))
  if(sum(!rownames(dat01) == samples$Basename) > 0) stop("Basenames don't match.")

  rownames(samples) = samples$Basename
  
  if(sum(!rownames(samples) == rownames(dat01)) > 0) stop("Error: rows don't match.")
  
  ## Average
  dat01 = aggregate(x = dat01, by = list(samples$SpeciesLatinName), FUN = function(x) return(mean(x, na.rm = T)))
  #dat01 = as.data.frame(dat01)
  colnames(dat01)[1] = "SpeciesLatinName"
  rownames(dat01) = dat01$SpeciesLatinName
  dat01 = dat01[, 2:ncol(dat01)]
  
  return(as.matrix(dat01))   
}



## Different titles control different outputs.
## To get averaged data from young samples (<3 years) only, use Young3_... This is one of the supplementary analysis

#mytitle = "Young3_Averaged_methCombined"
#mytitle = "Young5_NotMature_Averaged_methCombined"
mytitle = "Averaged_methCombined"
if(noExperiments == T) {
  mytitle = paste0("noExperiment_", mytitle)  
}

combinedsamples = combinedsamples[!is.na(combinedsamples$SpeciesLatinName), ]

if(grepl("Young|young", mytitle)) {    
  combinedsamples = combinedsamples[!is.na(combinedsamples$ConfidenceInAgeEstimate) & combinedsamples$ConfidenceInAgeEstimate >= 90 & !is.na(combinedsamples$Age),  ]
  
  agethreshold = sapply(strsplit(mytitle, "Young"), "[", 2)
  agethreshold = sapply(strsplit(agethreshold, "_"), "[", 1)
  if(agethreshold == "") {
    
    combinedsamples = combinedsamples[combinedsamples$Age < combinedsamples$averagedMaturity.yrs, ]
  } else {
    agethreshold = as.numeric(agethreshold)
    combinedsamples = combinedsamples[combinedsamples$Age < agethreshold, ]
  }
  
  if(grepl("NotMature", mytitle)) {
    combinedsamples = combinedsamples[combinedsamples$Age < combinedsamples$averagedMaturity.yrs, ] 
    
  }
  
  print(paste0("Young samples, oldest age is ", max(combinedsamples$Age)))
} 


## This should be the output from SeSAME normalization output, DNAm Beta values
load('./all_probes_all_samples_sesame.RData')
print(dim(dat0sesame))

## SESAME inverts x-y axes, so we transpose the matrix
dat0sesame = dat0sesame[!grepl("rs", dat0sesame$CGid), ]
dat0sesame = dat0sesame[, colnames(dat0sesame) %in% c("CGid", combinedsamples$Basename)]

cpgs = dat0sesame$CGid
dat0sesame = dat0sesame[, -1]
dat0sesame = t(dat0sesame)
colnames(dat0sesame) = cpgs


if(grepl("Young", mytitle)) {
  Averaged_methCombined = NoobAverageSpecies(combinedsamples, dat0sesame)
  #rm(samples, dat0sesame, NoobAverageSpecies)
  print(nrow(Averaged_methCombined))
  
  save(Averaged_methCombined, file = paste0("./", mytitle, ".RData"))
  
} else {
  Averaged_methCombined = NoobAverageSpecies(combinedsamples, NA, saveProof = "all")
  #rm(samples, dat0sesame, NoobAverageSpecies)
  print(nrow(Averaged_methCombined))
  
  
  ## Also save averaged data for specific tissue only, if needed.
  if(nrow(Averaged_methCombined) < 337) stop("Averaged_methCombined having less than 231 rows.")
  
  save(Averaged_methCombined, file = paste0("./", mytitle, ".RData"))
  
  Averaged_methCombined = NoobAverageSpecies(combinedsamples, NA, subsetVar = c("Tissue", "Liver"))
  save(Averaged_methCombined, file = paste0("./", "Liver_", mytitle, ".RData"))
  
  Averaged_methCombined = NoobAverageSpecies(combinedsamples, NA, subsetVar = c("Tissue", "Blood", "Spleen"), saveProof = "Blood")
  save(Averaged_methCombined, file = paste0("./", "Blood_", mytitle, ".RData"))
  
  Averaged_methCombined = NoobAverageSpecies(combinedsamples, NA, subsetVar = c("Tissue", "Skin", "Ear"))
  save(Averaged_methCombined, file = paste0("./", "Skin_", mytitle, ".RData"))
  
  Averaged_methCombined = NoobAverageSpecies(combinedsamples, NA, subsetVar = c("Tissue", "Muscle"))
  save(Averaged_methCombined, file = paste0("./", "Muscle_", mytitle, ".RData"))
  
  Averaged_methCombined = NoobAverageSpecies(combinedsamples, NA, subsetVar = c("Tissue", "Brain", "Astrocyte", "Cerebellum", "CortexEntorhinalis" ,"Cortex", "Hippocampus" , "Hypothalamus", "Striatum",  "SVZ", "Pituitary","nuclei_neuN_neg", "nuclei_neuN_pos"))
  save(Averaged_methCombined, file = paste0("./", "Brain_", mytitle, ".RData"))
  
}