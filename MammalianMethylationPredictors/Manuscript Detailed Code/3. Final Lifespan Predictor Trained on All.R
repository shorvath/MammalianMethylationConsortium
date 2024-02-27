# run master model in batch
rm(list=gc()); gc()


library(dplyr)
library(doParallel)
library(glmnet)


# load annotation data
anageCaesar <- read.csv(paste0(stuffcaesar, "anAgeUpdatedCaesarVersion49.csv"),
                        stringsAsFactors = FALSE)
source('./utilities_analysis.R')


#### Calibrations ###############
#################################

outcomeName = "logmaxAgeCaesar"
#outcomeName = "logaveragedMaturity.yrs"
#outcomeName = "logweightCaesar"
#outcomeName = "logGestation.Incubation..days."
pred.folds = 10
if(averagedCG == T) pred.folds = NULL
#selectedCGs = phyloResults$site[1:nSelect]

mytitle = "pred_AdjustLifespanAggressive_Overlap320K40K_DetectReady0.85_0.05FDR"

## in glmnet, 0.5 means it's Elastic Net
if(grepl("ridge|Ridge", mytitle)) {
    myalpha = 0  
} else {
    myalpha = 0.5
}

#################################################


### The 3 if statements, produce final lifespan predictor, species-aware predictor, and young speices only predictor, respectively
## load Averaged data produced by "Average CpGs by Species.R". Please run this file first to produce species-averaged data
## Since Mammalian lifespan predictor is trained on species-level data, we need to average all-sample data by species.
if(!grepl("SpeciesTissue", mytitle)) {
    load(paste0("./Averaged_methCombined.RData"))
    #Averaged_methCombined = Averaged_methCombined[!rownames(Averaged_methCombined) == "Balaena mysticetus", ]
    dat0 = Averaged_methCombined
    
    samples = data.frame(SpeciesLatinName = rownames(dat0))
    rm(Averaged_methCombined)
    
} else if(grepl("SpeciesTissue", mytitle)) {
    load(paste0("./SpeciesTissue_Averaged_methCombined.RData"))
    dat0 = Averaged_methCombined
    
    ## Filter out iPSC ES and Fibroblast
    dat0 = dat0[!grepl("iPSC|ES", rownames(dat0)), ]
    if(grepl("yumanensis|Yumanensis", mytitle)) {
        dat0 = dat0[!grepl("yumanensis", rownames(dat0)), ]
    }
    if(grepl("noHydroFibro", mytitle)) {
        dat0 = dat0[!grepl("Hydrochoerus hydrochaeris_Fibroblast", rownames(dat0)), ]
    }
    
    samples = data.frame(SpeciesTissue = rownames(dat0))
    samples$SpeciesTissue = as.character(samples$SpeciesTissue)
    samples$SpeciesLatinName = sapply(strsplit(samples$SpeciesTissue, "_"), '[', c(1))
}  else if(grepl("Young|young", mytitle)) {
    load(paste0("J./Young_Averaged_methCombined.RData"))
    dat0 = Averaged_methCombined
    samples = data.frame(SpeciesLatinName = rownames(dat0))
} 

## draw lifespan info from anageCaesar ###
samples = samples %>% left_join(anageCaesar, by = "SpeciesLatinName") 

## Adjust lifespan estimates based on Lu T. A. (2022) to be connsistent
samples$maxAgeCaesar = ifelse(samples$SpeciesLatinName %in% c("Homo sapiens", "Mus musculus"), 
                              samples$maxAgeCaesar, samples$maxAgeCaesar*1.3)

    
samples$logmaxAgeCaesar = log(samples$maxAgeCaesar)
samples$logweightCaesar = log(samples$weightCaesar)
samples$logaveragedMaturity.yrs = log(samples$averagedMaturity.yrs)
samples$logGestation.Incubation..days. = log(samples$Gestation.Incubation..days.)
samples$OrderFamily = paste0(samples$Order, "_", samples$Family)

## Shouldn't have NAs
if(sum(is.na(samples[, outcomeName])) > 0) {
    mystring = samples$SpeciesLatinName[is.na(samples[, outcomeName])]
    mystring = paste(mystring, collapse = " ")
    warning(paste0("NAs exist in outcome for ", mystring, ", filtered out."))
}
#samples = samples[!is.na(samples[, outcomeName]), ]
print(paste0("Total outcome NAs:", sum(is.na(samples[, outcomeName]))))
samples = samples[!is.na(samples[, outcomeName]), ]
##


if(grepl("Young|young", mytitle)) {
  dat0 = dat0[samples$SpeciesLatinName, ]
} else if(grepl("SpeciesTissue", mytitle)) {
  dat0 = dat0[samples$SpeciesTissue, ]
} else {
  dat0 = dat0[samples$SpeciesLatinName, ]
} 


################### Screening CGs ########################
## As described in Li. C. Z. et. al. (2024), filter CpGs by detection p-values
## Please use "CG Screening by DetectionP Value.R" to derive the results first
## The code in "CG Screening by DetectionP Value.R" takes in all Mammalian 40K Array detection p-value output from Sesame normalization pipeline
## And take the median p-value (FDR adjusted) of each probe per species, and keep those with 85% of species Median value >= 0.05
## Rationale: only want to keep probes that work in most of the species.
## In the end, load the "good" CpG list from the RDS file to be used for Elastic Net training code.
if(grepl("DetectReady", titleName)) {
    
    if(grepl("0.05FDR", mytitle)) {
        cpgs = readRDS(paste0(stuffcaesar, "JunoFork/DetectionP_median85_CGs_v11", myyuman, ".RDS"))
        print("Using FDR adjusted")
    }
     
    dat0 = dat0[, which(colnames(dat0) %in% c(cpgs))]
    print(paste0("Number of selected CGs: ", ncol(dat0), ". th=85"))
}

## No training-test set split. This is the final model.
samples$randos = 0
pred.folds = 1
produceCSV = T
speciesnumbers = data.frame(SpeciesLatinName = "LOO")

#########

gc()
#cores = detectCores()
#cl <- makeCluster((2), type = "PSOCK")
#registerDoParallel(cl)
#registerDoParallel(2)

print(paste0("dim(dat0)", dim(dat0)))
print(dim(samples))

train = as.data.frame(samples)    

## Train the model using glmnet
glmnet.Training.CV = cv.glmnet(as.matrix(dat0[!samples$randos == 1, ]), 
                               as.matrix(samples[!samples$randos == 1, outcomeName, drop = FALSE]), 
                               nfolds=10,
                               alpha=.5,family="gaussian")

glmnet.Training = glmnet(as.matrix(dat0[!samples$randos == 1, ]), 
                         as.matrix(samples[!samples$randos == 1, outcomeName, drop = FALSE]), 
                         family="gaussian", alpha=myalpha, nlambda=100)

# select the minimal lambda
lambda.min.glmnet.Training = glmnet.Training.CV$lambda.min
lambda.1se.glmnet.Training = glmnet.Training.CV$lambda.1se

# Model prediction
# train
#train = as.data.frame(samples)
#colnames(train)[1] = "Female"

if(produceCSV == T) {
 train$Y.pred <- as.numeric(predict(glmnet.Training, as.matrix(dat0), type="response",s=lambda.min.glmnet.Training))
 train$Y.pred.prob <- as.numeric(predict(glmnet.Training, as.matrix(dat0), type="response",s=lambda.min.glmnet.Training))
}

## Save coefficients of the final model
myfinal=data.frame(as.matrix(coef(glmnet.Training,s=lambda.min.glmnet.Training)))

colnames(myfinal) = speciesnumbers$SpeciesLatinName
saveRDS(myfinal, paste0("coefficients.RDS"))

    
if(T) {
    write.csv(train, paste0(titleName, ".csv"), row.names = FALSE)
}

rm(list = ls());gc()
