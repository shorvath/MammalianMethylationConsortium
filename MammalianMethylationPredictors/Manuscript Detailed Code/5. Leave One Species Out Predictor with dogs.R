
library(dplyr)
library(doParallel)
library(glmnet)

## Average dog methylation data

##########################################################################
############## Fit Dog Lifespan Predictor ONLY ############################ 
samples = read.csv("DogBreeds.csv",
                   stringsAsFactors = FALSE)

samples$logLifespanMedianHorvath = log(samples$LifespanMedianHorvath)
samples$logWeight.kg.avg = log(samples$Weight.kg.avg)

dat0 = readRDS("averaged_DogBreeds.RDS")

if(sum(!samples$DogBreed == rownames(dat0)) > 0) {
  stop("DogBreed must match dat0 rownames")
}

## Filter CpG
# This comes from filter CG R file
cpgs = readRDS(paste0("/DetectionP_median85_CGs_v11_AddYumanensis.RDS"))
dat0 = dat0[, colnames(dat0) %in% cpgs]
## Adjustment MULTIPLY by 1.3
samples$LifespanMedianHorvath = samples$LifespanMedianHorvath * 1.3

## Construct training test separation variable
outcomeName = "logLifespanMedianHorvath"
leavename = "DogBreed"
temp = data.frame(Order = unique(samples[, leavename])); colnames(temp) = leavename
temp$randos = 1:nrow(temp)
pred.folds = nrow(temp)
samples = samples %>% left_join(temp, by = leavename)


library(doParallel)
registerDoParallel(2)
print(paste0("dim(dat0)", dim(dat0)))
print(dim(samples))
test.set = foreach(i = 1:pred.folds, .combine = rbind, .packages = c('doParallel', 'glmnet')) %dopar% {
  
  manual.folds = NA 
  species.weights = NA
  
  # print(paste0("Weighting: ", species.weights[1]))
  lifespanpredictor = elasticNetTrain(dat0, samples, 
                                      outcomeName = outcomeName, selectLambda = "min", 
                                      leaveOneOut = i, leaveOutName = "randos", manual.folds = manual.folds, 
                                      species.weights = NA)
  
  toreturn = lifespanpredictor$test
  toreturn$Nselected = nrow(lifespanpredictor$glmnet.final.nonzero)
  return(toreturn)  
}

write.csv(test.set, paste0("./dog lifespan predictor.csv"), row.names = FALSE)



######################################################################################
############## Fit Combined Species + DogBreeds Predictor ############################ 
# leave one clade out

library(dplyr)
library(doParallel)
library(glmnet)
library(ggplot2)


SpeciesTissue = "AdjustLifespanAggressive_"
# SpeciesTissue = "NoMonotremata_"
# SpeciesTissue = "SpeciesTissue_"
# SpeciesTissue = "AddAge_"
# SpeciesTissue = ""

## This string controlls what the outcome is. In this case it's maximum lifespan
outcomeName = "logmaxAgeCaesar"
# outcomeName = "logmaxanAge"


## This string controlls what training-test split to use, leave one clade out?
## one species out? or leave out "small" mammals
## It's further explained in below ifelse statements, each correspond to a supplementary analysis
# leavename = "LeaveSmall"
# leavename = "Order"
# leavename = "OrderFamily"
leavename = "SpeciesLatinName"
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
  load('./Averaged_methCombined.RData')
}

dat0 = Averaged_methCombined[, !grepl("rs", colnames(Averaged_methCombined))] # filter out control probes

## Match KNN samples
dat0 = dat0[!rownames(dat0) %in% 
              c("Eulemur collaris", "Gerbillus cheesmani", "Sus scrofa minusculus", 
                "Cephalorhynchus hectori maui", "Rattus norvegicus domestica", "Galea leucoblepharum",
                "Paraechinus hypomelas", "Microgale drouhardi", "Microgale thomasi"), ]

if(SpeciesTissue == "NoMonotremata_") {
  temp = anageCaesar[anageCaesar$profiled & anageCaesar$Order == "Monotremata", "SpeciesLatinName"]
  dat0 = dat0[!rownames(dat0) %in%temp, ]
}

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


## Adjust lifespan estimates based on Lu T. A. (2022) to be connsistent
samples$maxAgeCaesar = ifelse(samples$SpeciesLatinName %in% c("Homo sapiens", "Mus musculus"), 
                              samples$maxAgeCaesar, samples$maxAgeCaesar*1.3)
samples$maxanAge = ifelse(samples$SpeciesLatinName %in% c("Homo sapiens", "Mus musculus"), 
                          samples$maxanAge, samples$maxanAge*1.3)


samples$logmaxAgeCaesar = log(samples$maxAgeCaesar)
samples$logweightCaesar = log(samples$weightCaesar)
samples$logaveragedMaturity.yrs = log(samples$averagedMaturity.yrs)
samples$logGestation.Incubation..days. = log(samples$Gestation.Incubation..days.)
samples$OrderFamily = paste0(samples$Order, "_", samples$Family)
samples$logmaxanAge = log(samples$maxanAge)


#### Add dog samples ######
dogs = read.csv("./DogBreeds.csv",
                stringsAsFactors = FALSE)

dogs$logLifespanMedianHorvath = log(dogs$LifespanMedianHorvath)
dogs$logWeight.kg.avg = log(dogs$Weight.kg.avg)
# ready to merge with dogs
dogs$SpeciesLatinName = dogs$DogBreed
dogs$maxAgeCaesar = dogs$LifespanMedianHorvath
dogs$logmaxAgeCaesar = dogs$logLifespanMedianHorvath
dogs$weightCaesar = dogs$Weight.kg.avg
dogs$logweightCaesar = dogs$logWeight.kg.avg
dogs = dogs[, c("SpeciesLatinName", "maxAgeCaesar", "logmaxAgeCaesar", "weightCaesar", "logweightCaesar")]

samples = samples[, c("SpeciesLatinName", "maxAgeCaesar", "logmaxAgeCaesar", "weightCaesar", "logweightCaesar")]
samples = rbind(samples, dogs)

dogdat0 = readRDS("./averaged_DogBreeds.RDS")

dogdat0 = dogdat0[, !grepl("rs", colnames(dogdat0))]
dogdat0 = dogdat0[, colnames(dogdat0) %in% cpgs]

if(sum(!samples$DogBreed == rownames(dogdat0)) > 0) {
  stop("DogBreed must match dat0 rownames")
}

dat0 = rbind(dat0, dogdat0)

#############################


if(outcomeName == "logmaxanAge") {
  dat0 = dat0[!is.na(samples$maxanAge), ]
  samples = samples[!is.na(samples$maxanAge), ]
  if(sum(!samples$SpeciesLatinName == rownames(dat0)) > 0) stop("LatinNames Don't Match")
}


titleName = paste0(mypath, "DogCombined_", outcomeName, "_", SpeciesTissue, 
                   leaveselect,"leaveOne", 
                   leavename, "_")
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


