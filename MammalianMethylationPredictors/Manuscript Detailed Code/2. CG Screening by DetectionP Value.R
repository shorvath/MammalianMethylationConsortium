


## As described in Li. C. Z. et. al. (2024), filter CpGs by detection p-values
## This code takes in all Mammalian 40K Array detection p-value output from Sesame normalization pipeline
## And take the median p-value (FDR adjusted) of each probe per species, and keep those with 85% of species Median value >= 0.05
## Rationale: only want to keep probes that work in most of the species.
## In the end, save the file to an RDS file to be used for Elastic Net training code.

  
load("./detectionP_combined_sesame.RData")
detectionP_combined_sesame = detectionP_combined_sesame[!grepl("rs", detectionP_combined_sesame$CGid), ]

if(grepl("EPIC450", titleName)) {
  epic = read.csv("~/StuffCaesar/utilities/ProbesSharedWithEPIC.CanonicalManifestMinfi_overlapwith450k.csv")
  detectionP_combined_sesame = detectionP_combined_sesame[detectionP_combined_sesame$CGid %in% epic$IlmnID, ]
} else if(grepl("EPIC", titleName))  {
  epic = read.csv("~/StuffCaesar/utilities/ProbesSharedWithEPIC.CanonicalManifestMinfi.csv")
  detectionP_combined_sesame = detectionP_combined_sesame[detectionP_combined_sesame$CGid %in% epic$IlmnID, ]
}

#cpgs = detectionP_combined_sesame$CGid
#detectionP_combined_sesame = detectionP_combined_sesame[, -1]
#detectionP_combined_sesame = t(detectionP_combined_sesame)
#colnames(detectionP_combined_sesame) = cpgs
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
combinedsamples = combinedsamples[!grepl("N71\\.|N79\\.|N80\\.|N82\\.|N83\\.|N85\\.|N128", combinedsamples$Folder), ]

## ADD only N141
combinedsamples = combinedsamples[!grepl("N138\\.|N139\\.|N140\\.|N142\\.", combinedsamples$Folder), ]

combinedsamples = combinedsamples[!grepl("iPSC|ES", combinedsamples$Tissue), ]
combinedsamples = combinedsamples[combinedsamples$CanBeUsedForAgingStudies == "yes", ]
combinedsamples = combinedsamples[!grepl("Tt x Tt gilli", combinedsamples$SpeciesLatinName), ]

if(grepl("Yumanensis", mytitle)) {
  combinedsamples = combinedsamples[!grepl("yumanensis", combinedsamples$SpeciesLatinName), ]
}
if(grepl("noHydroFibro", mytitle)) {
  combinedsamples = combinedsamples[!(combinedsamples$SpeciesLatinName == "Hydrochoerus hydrochaeris" & 
                                        combinedsamples$Tissue == "Fibroblast"), ]
}

if(grepl("noExperiment", mytitle)) {
  combinedsamples = read.csv("~/StuffCaesar/utilities/combinedsamples_noExperiments.csv", stringsAsFactors = F)
}

###
#combinedsamples = combinedsamples[combinedsamples$Basename %in% colnames(detectionP_combined_sesame), ]
###

xy <- prepDat0(combinedsamples, detectionP_combined_sesame, filterCanBeUsed = TRUE, 
               filterAgeConfidence = FALSE, filterCharacteristics = NULL, filterFolders = NULL,
               filterFewSpecies = FALSE)
combinedsamples <- xy$ys
detectionP_combined_sesame <- xy$x
rm(xy); gc()
detectionP_combined_sesame = 
  detectionP_combined_sesame[, which(colnames(detectionP_combined_sesame) %in% colnames(dat0))]

combinedsamples = combinedsamples %>% left_join(anageCaesar, by = "SpeciesLatinName")

#if(grepl("DetectionP", titleName)) {
#    combinedsamples$SpeciesLatinName = paste0(combinedsamples$SpeciesLatinName, "_", combinedsamples$Tissue)
#}

detectionP_combined_sesame = as.data.frame(detectionP_combined_sesame)
cpgs = colnames(detectionP_combined_sesame)
detectionP_combined_sesame$SpeciesLatinName = combinedsamples$SpeciesLatinName

detectionP_combined_sesame = detectionP_combined_sesame %>% 
  group_by(SpeciesLatinName) %>%
  summarise_at(cpgs, median)
detectionP_combined_sesame = as.data.frame(detectionP_combined_sesame)
rownames(detectionP_combined_sesame) = detectionP_combined_sesame$SpeciesLatinName
detectionP_combined_sesame = detectionP_combined_sesame[, -which(colnames(detectionP_combined_sesame) == "SpeciesLatinName")]
detectionP_combined_sesame = as.data.frame(detectionP_combined_sesame)

cpgs = colnames(detectionP_combined_sesame)
basenames = rownames(detectionP_combined_sesame)
detectionP_combined_sesame = as.matrix(detectionP_combined_sesame)
temp = p.adjust(detectionP_combined_sesame, method = "fdr")
detectionP_combined_sesame = matrix(data = temp, nrow = nrow(detectionP_combined_sesame), byrow = FALSE)
colnames(detectionP_combined_sesame) = cpgs
rownames(detectionP_combined_sesame) = basenames

filterrows = ceiling(nrow(detectionP_combined_sesame) * 0.85)
#filterrows = 140
cpgs = colnames(detectionP_combined_sesame)
cpgs = cpgs[apply(detectionP_combined_sesame, 2, function(x) return(sum(x < 0.05) >= filterrows))]

# FILTER OUT SNPs
cpgs = cpgs[!grepl("rs", cpgs)]

saveRDS(cpgs, file = paste0("YOUR_PATH/DetectionP_median85_CGs_v11.RDS"))
  
