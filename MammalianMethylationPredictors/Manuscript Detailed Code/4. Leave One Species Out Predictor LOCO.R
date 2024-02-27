
######################################################################################
############## Leave One Clade Out Predictor ############################ 
# leave one clade out

library(dplyr)
library(doParallel)
library(glmnet)


## This string controlls what the outcome is. In this case it's maximum lifespan
outcomeName = "logmaxAgeCaesar"
# outcomeName = "logmaxanAge"


## This string controlls what training-test split to use, leave one clade out?
## one species out? or leave out "small" mammals
## It's further explained in below ifelse statements, each correspond to a supplementary analysis
# leavename = "LeaveSmall"
leavename = "Order"
# leavename = "OrderFamily"
# leavename = "SpeciesLatinName"
# leavename = "LeaveCetacea"
leaveselect = 0

source("./utilities_analysis.R")

# anAge data
anageCaesar <- read.csv("anAgeUpdatedCaesarVersion49.csv",
                        stringsAsFactors = FALSE)

## load Averaged data produced by "Average CpGs by Species.R". Please run this file first to produce species-averaged data
## Since Mammalian lifespan predictor is trained on species-level data, we need to average all-sample data by species.
if(SpeciesTissue == "SpeciesTissue_") {
  load('./SpeciesTissue_Averaged_methCombined.RData')
} else {
  load('./AveragedbySpecies/Version8_March2022/Averaged_methCombined.RData')
}

dat0 = Averaged_methCombined[, !grepl("rs", colnames(Averaged_methCombined))]

#Averaged_methCombined = Averaged_methCombined[!rownames(Averaged_methCombined) == "Balaena mysticetus", ]

## Match KNN samples
## Note that for KNN, the TimeTree project lacks a few species. 
## Therefore we remove these species in our training set as well to have a fair comparison
dat0 = dat0[!rownames(dat0) %in% 
              c("Eulemur collaris", "Gerbillus cheesmani", "Sus scrofa minusculus", 
                "Cephalorhynchus hectori maui", "Rattus norvegicus domestica", "Galea leucoblepharum",
                "Paraechinus hypomelas", "Microgale drouhardi", "Microgale thomasi"), ]


## As described in Li. C. Z. et. al. (2024), filter CpGs by detection p-values
## Please use "CG Screening by DetectionP Value.R" to derive the results first
## The code in "CG Screening by DetectionP Value.R" takes in all Mammalian 40K Array detection p-value output from Sesame normalization pipeline
## And take the median p-value (FDR adjusted) of each probe per species, and keep those with 85% of species Median value >= 0.05
## Rationale: only want to keep probes that work in most of the species.
## In the end, load the "good" CpG list from the RDS file to be used for Elastic Net training code.
cpgs = readRDS(paste0("./DetectionP_median85_CGs_v11_AddYumanensis.RDS"))

dat0 = dat0[, colnames(dat0) %in% cpgs]

## Matching genotype data to phenotype annotations
samples = data.frame(SpeciesLatinName = rownames(dat0))
rm(Averaged_methCombined)
samples = samples %>% left_join(anageCaesar, by = "SpeciesLatinName") 

  
## Calibration to be consistent with A.T. Lu et al. (2019)
samples$maxAgeCaesar = ifelse(samples$SpeciesLatinName %in% c("Homo sapiens", "Mus musculus"), 
                              samples$maxAgeCaesar, samples$maxAgeCaesar*1.3)
samples$maxanAge = ifelse(samples$SpeciesLatinName %in% c("Homo sapiens", "Mus musculus"), 
                          samples$maxanAge, samples$maxanAge*1.3)
##

samples$logmaxAgeCaesar = log(samples$maxAgeCaesar)
samples$logweightCaesar = log(samples$weightCaesar)
samples$logaveragedMaturity.yrs = log(samples$averagedMaturity.yrs)
samples$logGestation.Incubation..days. = log(samples$Gestation.Incubation..days.)
samples$OrderFamily = paste0(samples$Order, "_", samples$Family)
samples$logmaxanAge = log(samples$maxanAge)


if(outcomeName == "logmaxanAge") {
  dat0 = dat0[!is.na(samples$maxanAge), ]
  samples = samples[!is.na(samples$maxanAge), ]
  if(sum(!samples$SpeciesLatinName == rownames(dat0)) > 0) stop("LatinNames Don't Match")
}

## Creating a new variable to split training-test set folds
temp = data.frame(Order = unique(samples[, leavename])); colnames(temp) = leavename
temp$randos = 1:nrow(temp)
pred.folds = nrow(temp)
samples = samples %>% left_join(temp, by = leavename)


## In the paper, for LOCO, leaveselect is set to 2, keep only 2 species from large Mammalian Orders in training set
if(leaveselect > 0 & leavename == "Order") {
  # set.seed(123)
  if(leaveselect == 2) {
    for (i in 1:length(unique(samples[, "randos"]))) {
      # myindices = which(samples[, "randos"] == i)
      myindices = samples[samples[, "randos"] == i, ]
      if(nrow(myindices) >= 20) {
        indextemp = which(myindices[, outcomeName] == max(myindices[, outcomeName]) | 
                            myindices[, outcomeName] == min(myindices[, outcomeName]))
        # indextemp = c(indextemp, which.min(abs(myindices[, "logmaxAgeCaesar"] - median(myindices[, "logmaxAgeCaesar"]))))
        
        myindices = myindices[indextemp, , drop = FALSE]
        # myindices = sample(myindices, leaveselect, replace = FALSE)
        samples$randos[samples$SpeciesLatinName %in% myindices$SpeciesLatinName] = 0
      }
    }
  } 
  
## Leave small species out to produce Supplementary S14
} else if(grepl("LeaveSmall", leavename)) {
  weightThreshold = 150
  samples$randos = ifelse(samples$weightCaesar < weightThreshold, 1, 2)
  myfdr = paste0("weightCut", weightThreshold)
  pred.folds =  1
} else if(grepl("Leave", leavename)) {
  tempname = sapply(strsplit(leavename, "Leave"), "[", c(2))
  samples$randos = NA
  samples$randos = ifelse(samples$Order == tempname, 1, 2)
  pred.folds =  1
} 


rm(temp)


## Using a parallelized code to speed up the training process.
## Connvenient function, elasticNetTrain() written to train a glmnet object in accordance to the genotype
## data, outcome, and splits
## For details please see sourced R code source("./utilities_analysis.R")

registerDoParallel(2)
print(paste0("dim(dat0)", dim(dat0)))
print(dim(samples))
test.set = foreach(i = 1:pred.folds, .combine = rbind, .packages = c('doParallel', 'glmnet')) %dopar% {
  
  manual.folds = NA 
  
  species.weights = NA
  
  # This function uses glmnet net to train Elastic Net models 
  lifespanpredictor = elasticNetTrain(dat0, samples, 
                                      outcomeName = outcomeName, selectLambda = "min", 
                                      leaveOneOut = i, leaveOutName = "randos", manual.folds = manual.folds, 
                                      species.weights = NA)
  
  toreturn = lifespanpredictor$test
  toreturn$Nselected = nrow(lifespanpredictor$glmnet.final.nonzero)
  return(toreturn)  
}


write.csv(test.set, paste0("PATH to Results"), row.names = FALSE)


