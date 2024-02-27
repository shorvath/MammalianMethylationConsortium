
################################################################################
#    Transform dat0 methylation data and filter                                #
#                                                                              #
################################################################################
prepDat0 <- function(samples = NULL, normalized = NULL, filterCanBeUsed = TRUE, filterAgeConfidence = TRUE, 
                     filterCharacteristics = NULL, filterFolders = NULL, filterAge = NULL, filterFewSpecies = NULL,
                     subsetVariable = c(NA, NA), filterSpeciesTissue = NULL, filterFew = NULL) {
  
  library(dplyr)
  if(!is.null(normalized)) {
    cpgs = normalized$CGid
    dat0 = normalized[, -1]
    dat0 = t(dat0)
    colnames(dat0) = cpgs
  }
  
  if(class(samples$Female) == "character") 
    samples$Female = ifelse(samples$Female == "1", 
                            1, ifelse(samples$Female == "0", 0, NA))
  
  rownames(samples) = samples$Basename
  
  if(!is.null(normalized)) {
    if(!nrow(samples) == nrow(dat0)) print("Warning: Rows don't match.")
    if(sum(!samples$Basename %in% rownames(dat0)) > 0) print("Warning: some basenames don't exits in dat0.")
    #samples = samples[rownames(dat0), ]
  }
  
  if(filterCanBeUsed == TRUE) {
    samples = samples[!is.na(samples$CanBeUsedForAgingStudies), ]
    samples = samples[samples$CanBeUsedForAgingStudies == "yes", ]
  }
  
  if(!is.null(filterAgeConfidence) & !is.na(filterAgeConfidence) & !filterAgeConfidence==F) {
    
    if(!is.numeric(filterAgeConfidence)) {filterAgeConfidence = 90}
    
    print(paste0("Filtering AgeConfidence by ", filterAgeConfidence, "%"))
    samples = samples[!is.na(samples$Age), ]
    samples = samples[samples$ConfidenceInAgeEstimate >= filterAgeConfidence & !is.na(samples$ConfidenceInAgeEstimate), ]
    #dat0 = dat0[rownames(samples), ]
  }
  
  if(!is.null(filterCharacteristics)) {
    for (i in 1:length(filterCharacteristics)) {
      samples = samples[!is.na(samples[, filterCharacteristics[i]]), ]
    }
  }
  
  # filterFolder can be a vector of folders
  if(!is.null(filterFolders)) {
    samples = samples[samples$Folder %in% filterFolders, ]
  }
  
  # filter age
  if(!is.null(filterAge)) {
    samples = samples[samples$Age < filterAge[2] & samples$Age > filterAge[1], ]
  }
  
  # filter too few Species samples for methylation rate analysis
  if(!is.null(filterFewSpecies)) {
    keptSpecies <- samples %>% group_by(SpeciesLatinName) %>% summarise(n = n()) %>% filter(n >= filterFewSpecies)
    keptSpecies <- keptSpecies$SpeciesLatinName
    samples <- samples[samples$SpeciesLatinName %in% keptSpecies, ]
  }
  if(!is.null(filterSpeciesTissue)) {
    samples$SpeciesTissue = paste0(samples$SpeciesLatinName, "_", samples$Tissue)
    keptSpecies <- samples %>% group_by(SpeciesTissue) %>% summarise(n = n()) %>% filter(n >= filterSpeciesTissue)
    keptSpecies <- keptSpecies$SpeciesTissue
    samples <- samples[samples$SpeciesTissue %in% keptSpecies, ]
  }
  if(!is.null(filterFew)) {
    samples$myvar = samples[, filterFew[1]]
    keptSpecies <- samples %>% group_by(myvar) %>% summarise(n = n()) %>% filter(n >= as.numeric(filterFew[2]))
    keptSpecies <- keptSpecies$myvar
    samples <- samples[samples$myvar %in% keptSpecies, ]
    samples = samples[, -which(colnames(samples) == "myvar")]
  }
  
  # Subset on a variable (e.g. Tissues)
  if(!is.na(subsetVariable[1])) {
    samples = samples[!is.na(samples[, subsetVariable[1]]), ]
    samples = samples[samples[, subsetVariable[1]] %in% subsetVariable, ]
  }
  
  # filter methylation data accordingly to match the samplesheet
  if(!is.null(normalized)) {
    dat0 = dat0[rownames(samples), ]
  } else {
    dat0 = NULL
  }
  
  toreturn <- structure(list(x = dat0,
                             ys = samples),
                        class = "cpgReady")
  
  return(toreturn)
}

################################################################################
#    Epigenetic Clock Plots                                                    #
#                                                                              #
################################################################################
epiClockPlot <- function(glmnet.Training.CV = NULL, lambda.glmnet.Training = NULL, train = NULL, test = NULL, myobject = NULL,
                         addColors = NULL, outcomeName = NULL) {
  library(dplyr)
  library(WGCNA)
  if(!is.null(myobject)) {
    glmnet.Training.CV = myobject$glmnet.Training.CV 
    glmnet.Training = myobject$glmnet.Training
    lambda.glmnet.Training = myobject$lambda.glmnet.Training
    train = myobject$train
    test = myobject$test
  }
  if(is.null(outcomeName)) outcomeName = colnames(train)[2]
  
  # if(!"Age" %in% colnames(train) | !"Age" %in% colnames(test)) {
  #   train$Age = train$relativeAge
  #   test$Age = test$relativeAge
  # }
  
  if(!is.null(addColors)) {
    if(addColors[1] == TRUE) {
      temp = test %>% group_by(SpeciesLatinName) %>% summarize(n = n()) %>% arrange(n)
      if(nrow(temp) > 5) temp = temp$SpeciesLatinName[1:5]
      else temp = temp$SpeciesLatinName
    } else {
      speciesColors.test = addColors
      colorfactor2 = as.factor(test$DogBreed)
      temp = levels(colorfactor2)
      todisplay2 = match(temp, levels(colorfactor2))
    }
    #speciesColors.train = addColors[rownames(train), "SpeciesLatinName"]
    #colorfactor1 = as.factor(train$SpeciesLatinName)
    #todisplay1 = match(temp, levels(colorfactor1))
    #if(!length(todisplay1) == 3) stop("Color failed")
    #speciesColors.test = addColors[rownames(test), "SpeciesLatinName"]
    #colorfactor2 = as.factor(test$SpeciesLatinName)
    #todisplay2 = match(temp, levels(colorfactor2))
  } 
  # Plot function
  #par(mfrow=c(2,2))
  if(!is.null(addColors)) {
    layout.matrix = matrix(c(1,2,3,4, 5, 5), ncol=2, byrow=TRUE)
    layout.heights = c(4, 4, 1.5)
  } else {
    layout.matrix = matrix(c(1,2,3,4), ncol=2, byrow=TRUE)
    layout.heights = c(4, 4)
  }
  layout(layout.matrix, heights=layout.heights )
  par(mai=c(1, .8, .5, .5))
  plot(glmnet.Training.CV, main=paste0('optimal lambda=', round(lambda.glmnet.Training, digits=3)))
  title('A',adj=0,font=2,cex.main=2)
  par(mai=c(1, .8, .5, .5))
  plot(glmnet.Training, label = TRUE)
  title('B',adj=0,font=2,cex.main=2)
  
  # Function for dealing with undefined correlations
  Corr <<- function(a, b, use = "p") {
    mycor = ifelse(is.na(cor(a, b, use = use)), 0, cor(a, b, use = use))
    return(mycor)
  }
  # Model fit
  if(!is.null(addColors)) {
    par(mai=c(.8, .8, .5, .5))
    main.txt=paste0('Training Set, ')
    verboseScatterplot(train$Y.pred, train[, outcomeName], main = main.txt, col = 1,
         xlab = paste0("Predicted ", outcomeName), ylab = outcomeName, corFnc = "Corr")
    # legend(x = c(0, 25), y = c(55, 85),  y.intersp = 0.45, x.intersp = .4,
    #        legend = c("Mus musculus", "Homo sapiens", "Canis lupus familiaris"),
    #        pch = 1,  title="Top Occurring Species",
    #        col = todisplay1, horiz=FALSE, cex=1)
    abline(0,1)
    title('C', adj=0,font=2,cex.main=2)
    #
    par(mai=c(.8, .8, .5, .5))
    main.txt=paste0('Test Set, ')
    verboseScatterplot(test$Y.pred, test[, outcomeName], main = main.txt, col = colorfactor2,
                       xlab = paste0("Predicted ", outcomeName), ylab = outcomeName, corFnc = "Corr")
    # text(test$Y.pred[test$Age / test$Y.pred > 2], 
    #      test$Age[test$Age / test$Y.pred > 2],  
    #      test$Basename[test$Age / test$Y.pred > 2], cex=0.6, pos=1, col="red")
    # legend(x = c(0.1, .4), y = c(.45, .7), inset=.02, y.intersp = 0.45, x.intersp = .4,
    #        legend = temp,
    #        pch = 1, title="Top Occurring Species",
    #        col = todisplay2, horiz=FALSE, cex=1)
    # legend("topleft", inset=.02,
    #        legend = temp,
    #        pch = 1, title="Top Occurring Species",
    #        col = todisplay2, text.width=rep(0.1, length(temp)))
    abline(0,1)
    title('D',adj=0,font=2,cex.main=2)
    #
    par(mai=c(0,0,0,0))
    plot.new()
    legend(x = "center",legend = temp, x.intersp = 0.3,
           pch = 1, title="Top Occurring Species in Test Set", cex = 1.2,
           col = todisplay2, horiz = TRUE, text.width=rep(1.2, length(temp)))
  } else {
    par(mai=c(1, .8, .5, .5))
    main.txt=paste0('Training Set, ')
    verboseScatterplot(train$Y.pred, train[, outcomeName], main = main.txt,
                       xlab = paste0("Predicted ", outcomeName), ylab = outcomeName, corFnc = "Corr")
    abline(0,1)
    title('C', adj=0,font=2,cex.main=2)
    #
    par(mai=c(1, .8, .5, .5))
    main.txt=paste0('Test Set, ')
    verboseScatterplot(test$Y.pred, test[, outcomeName], main = main.txt,
                       xlab = paste0("Predicted ", outcomeName), ylab = outcomeName, corFnc = "Corr")
    # text(test$Y.pred[test$Age / test$Y.pred > 2], 
    #      test$Age[test$Age / test$Y.pred > 2],  
    #      test$Basename[test$Age / test$Y.pred > 2], cex=0.6, pos=1, col="red")
    abline(0,1)
    title('D',adj=0,font=2,cex.main=2)
    #
  }
}

################################################################################
#    Epigenetic Clock Analysis                                                 #
#                                                                              #
################################################################################
elasticNetTrain <- function(x, y, outcomeName = "Age", trainPROPORTION = .7, 
                            SEED = NULL, NFOLD=10,  ALPHA=0.5, selectLambda = "min",
                            leaveOneOut = NULL, leaveOutName = "SpeciesLatinName", 
                            manual.folds = NA, species.weights = NA, family = "gaussian",
                            filterNAs = FALSE) {
  if(family == "gaussian") {
    mytype = "response"
  } else {
    mytype = "class"
  }
  if(is.null(leaveOneOut)) {
    # set RNG type to be the same, in case an older R version is used
    # NOTE: R 3.6.0 USES A DIFFERENT RNG!
    if(!is.null(SEED))
      set.seed(SEED, kind = "Mersenne-Twister")
    trains <- sample(nrow(x), ceiling(trainPROPORTION*nrow(x)), replace = FALSE)
    # Removing random seed, good habbit. Random seed can be a dangerous thing
    rm(.Random.seed, envir=globalenv())
  } else {
    trains <- 1:nrow(x)
    if(filterNAs == TRUE) {
      trains <- trains[!y[, leaveOutName] %in% leaveOneOut & !is.na(y[, outcomeName])]
    } else {
      trains <- trains[!y[, leaveOutName] %in% leaveOneOut]
    }
  }

  # train.x <- x[trains,]
  # train.y <- y[trains, , drop = FALSE]
  # test.x <- x[-trains, ]
  # test.y <- y[-trains, , drop = FALSE]

  library(glmnet)
  library(dplyr)
  # Makevars
  # if(!is.na(manual.folds)) {
  #   manual.folds = y[trains, "SpeciesLatinName", drop = FALSE] %>% left_join(manual.folds, by = "SpeciesLatinName")
  #   NFOLD = NA
  #   manual.folds = manual.folds$folds
  # } 
  
  if(!is.na(manual.folds[1])) {
    manual.folds = y[trains, "foldid"]
  }
  
  if(typeof(species.weights) == "logical") {
    if(!is.na(species.weights))
      species.weights = y[trains, "species.weights", drop = TRUE]
  } 
  # else if(!is.na(species.weights)) {
  #   temp = y[trains, "SpeciesLatinName", drop = FALSE]
  #   temp = temp %>% left_join(species.weights, by = "SpeciesLatinName")
  #   species.weights = temp$n; rm(temp)
  # }
 
  #ALPHA=1#lasso
  #ALPHA=0#ridge
  #set.seed(SEED + 999, kind = "Mersenne-Twister")
  if(!is.na(manual.folds)) {
    glmnet.Training.CV = cv.glmnet(as.matrix(x[trains,]), as.matrix(y[trains, outcomeName, drop = FALSE]), 
                                   family = family, alpha = ALPHA, nfolds = NFOLD, foldid = manual.folds)
  } else if(sum(is.na(species.weights)) == 0) {
    glmnet.Training.CV = cv.glmnet(as.matrix(x[trains,]), as.matrix(y[trains, outcomeName, drop = FALSE]), 
                                   family = family, alpha = ALPHA, nfolds = NFOLD, weights = species.weights)
  } else {
    glmnet.Training.CV = cv.glmnet(as.matrix(x[trains,]), as.matrix(y[trains, outcomeName, drop = FALSE]), 
                                   family = family, alpha = ALPHA, nfolds = NFOLD)
  }
  #rm(.Random.seed, envir=globalenv())
  
  # Select Lambda
  if(selectLambda == "min") {
    lambda.glmnet.Training = glmnet.Training.CV$lambda.min
  } else {
    lambda.glmnet.Training = glmnet.Training.CV$lambda.1se
  }
  glmnet.Training = glmnet(as.matrix(x[trains,]), as.matrix(y[trains, outcomeName, drop = FALSE]), 
                           family = family, alpha = ALPHA, nlambda=100)
  
  # Predict Validation Set
  
  # Train set
  train = data.frame(y[trains, , drop = FALSE], is.train=1)
  train$Y.pred=as.numeric(predict(glmnet.Training, as.matrix(x[trains, , drop = FALSE]), type=mytype, s=lambda.glmnet.Training))
  #train$Y.pred.prob <- as.numeric(predict(glmnet.Training, as.matrix(train.x), type="response",s=lambda.min.glmnet.Training))
  
  # Test set
  test=data.frame(y[-trains, , drop = FALSE], is.train=0)
  test$Y.pred=as.numeric(predict(glmnet.Training, as.matrix(x[-trains, , drop = FALSE]), type=mytype,s=lambda.glmnet.Training))
  if(family == "binomial")
    test$Y.pred.prob <- as.numeric(predict(glmnet.Training,  as.matrix(x[-trains, , drop = FALSE]), type="response",s=lambda.glmnet.Training))
  
  glmnet.final=data.frame(as.matrix(coef(glmnet.Training,s=lambda.glmnet.Training)))
  names(glmnet.final)='beta'
  glmnet.final$var=rownames(glmnet.final)
  #glmnet.final=subset(glmnet.final,select=c(var,beta))
  glmnet.final.nonzero=subset(glmnet.final,abs(beta)>0)
  
  toreturn <- structure(list(glmnet.Training.CV = glmnet.Training.CV,
                             glmnet.Training = glmnet.Training,
                             lambda.min = glmnet.Training.CV$lambda.min,
                             lambda.1se = glmnet.Training.CV$lambda.1se,
                             train = train,
                             test = test,
                             glmnet.final.nonzero = glmnet.final.nonzero),
                        class = "elasticnet")
  
  return(toreturn)
}

################################################################################
#  Seperate SpeciesLatinName strings into ToL format  (Genus_species)          #
#                                                                              #
################################################################################
as.labelsCaesar <- function(SpeciesLatinName, noWarnings = T) {
  
  if(length(SpeciesLatinName) == 0) {return("")}
  
  if(is.factor(SpeciesLatinName)) SpeciesLatinName = as.character(SpeciesLatinName)
  j = 1
  totalLength = length(SpeciesLatinName)
  for (i in 1:totalLength) {
    temp = unlist(strsplit(SpeciesLatinName[i], "[^a-zA-Z]"))
    if(length(temp) > 2 & noWarnings == F) {
      if(j == 1)  {
        cat(paste0("WARNING: species name \"", SpeciesLatinName[i], "\","))
      } else {
        cat(paste0("\"", SpeciesLatinName[i], "\","))
      }
      j = j + 1
    }
    if(i == totalLength) cat("might have included a subspecies name. \n")
  }
  
  temp = sapply(strsplit(SpeciesLatinName, "[^a-zA-Z]"), '[', c(1,2))
  temp = t(temp) 
  
  temp[,1] = paste0(toupper(substr(temp[,1], 1, 1)), substr(temp[,1], 2, nchar(temp[,1])))
  temp[,2] = paste0(tolower(substr(temp[,2], 1, 1)), substr(temp[,2], 2, nchar(temp[,2])))
  
  labelsCaesar = paste(temp[,1], temp[,2], sep = "_")
  
  return(labelsCaesar)
}

################################################################################
#  Trim the Tree to have only nodes existing in the data                       #
#                                                                              #
################################################################################

trimPhyloTree <- function(samples = NULL, tree = NULL) {
  
  library(ape)
  # Create necessary tree lengths and species labels if not existing already
  if(is.null(tree$edge.length)) {
    tree = compute.brlen(tree, method = 'Grafen')
  }
  
  if(! "labelsCaesar" %in% colnames(samples) & ! "SpeciesLatinName" %in% colnames(samples)) {
    samples$SpeciesLatinName = rownames(samples)
  }
  if(! "labelsCaesar" %in% colnames(samples)) {
    samples$labelsCaesar <- as.labelsCaesar(samples$SpeciesLatinName)
  }
  
  rownames(samples) = samples$labelsCaesar
  dat.dropped = samples[samples$labelsCaesar %in% tree$tip.label, ]
  rownames(dat.dropped) = dat.dropped$labelsCaesar
  tree.dropped = keep.tip(tree, dat.dropped$labelsCaesar)
  dat.dropped = dat.dropped[tree.dropped$tip.label, ]
  return(list(dat.dropped, tree.dropped))
  
}

################################################################################
#  Correct SpeciesLatinName Typos for Tree                                      #
#                                                                              #
################################################################################

correct.names <- function(samples) {
  
  ## Naming conventions
  # Ake updates
  # Aonyx cinereus
  samples$SpeciesLatinName[samples$SpeciesLatinName == "Amblonyx cinereus"] = "Aonyx cinereus"
  # Ceratotherium simum simum
  samples$SpeciesLatinName[samples$SpeciesLatinName == "Equus ferus caballus"] = "Equus caballus"
  # Gazella granti
  samples$SpeciesLatinName[samples$SpeciesLatinName == "Gazella granti"] = "Nanger granti"
  # Marmota flaviventer
  samples$SpeciesLatinName[samples$SpeciesLatinName == "Marmota flaviventer"] = "Marmota flaviventris"
  # Marmota flaviventer
  samples$SpeciesLatinName[samples$SpeciesLatinName == "Macropus ruforiseus"] = "Macropus rufogriseus"
  # Marmota flaviventer
  samples$SpeciesLatinName[samples$SpeciesLatinName == "Peromycus polionotus"] = "Peromyscus polionotus"
  # Josh Names Corrected 12/14/2020
  # Source: word doc in StillMissing, Damaraland mole-rat.
  samples$SpeciesLatinName[samples$SpeciesLatinName == "Erethizon dorsatum"] = "Erethizon dorsatus"
  samples$SpeciesLatinName[samples$SpeciesLatinName == "Fukomys damarensis"] = "Cryptomys damarensis"
  samples$SpeciesLatinName[samples$SpeciesLatinName == "Nannospalax ehrenbergi"] = "Spalax ehrenbergi"
  samples$SpeciesLatinName[samples$SpeciesLatinName == "Mesoplon bidens"] = "Mesoplodon bidens"
  
  return(samples)
}

################################################################################
#  Correct some DogBreed name discrepancies                                    #
#  From finalFile names to -> Pedigree file                                    #
################################################################################

dogname_synchronizer <- function(namestring) {
  namestring[namestring == "Cocker Spaniel"] = "American Cocker Spaniel"
  namestring[namestring == "Flat-Coated Retriever"] = "Flat-coated Retriever"
  
  namestring[namestring == "Mastiff"] = "English Mastiff"
  namestring[namestring == "Miniature Poodle"] = "Poodle - Miniature"
  namestring[namestring == "Pug"] = "Pug Dog"
  namestring[namestring == "Standard Poodle"] = "Poodle - Standard"
  namestring[namestring == "Toy Poodle"] = "Poodle - Toy"
  
  return(namestring)
}

################################################################################
#  rmcorr function                                                             #
#  From finalFile names to -> Pedigree file                                    #
################################################################################
F1_ewas=function(CG=NULL, df=NULL, outcomeName = "Age"){
  library(rmcorr)
  df$Y = df[, outcomeName]
  df=df[, c("Basename", "Y", "SpeciesLatinName")]
  df$SpeciesLatinName1=as.numeric(as.factor(df$SpeciesLatinName))
  df$SpeciesLatinName1=as.factor(df$SpeciesLatinName1)
  
  df$CpG=CG
  m0=rmcorr(SpeciesLatinName1, Y, CpG, dataset=df)
  #m0=rmcorr(df$SpeciesLatinName1, df[, outcomeName], df$CpG)
  m.corr=cor.test(df$Y, df$CpG)
  x.out=data.frame(rmcor=m0$r,  p.rmcor=m0$p,  df=m0$df)
  return(x.out)
}
#

################################################################################
#  Correct Mappability Names                                                   #
#  From stuck-together names -> SpeciesLatinName                               #
################################################################################
unStuckSpeciesLatinNames <- function(SpeciesNameInput) {
  temp = gsub('([a-z])(?=[A-Z])','\\1 ', SpeciesNameInput, perl=T)
  return(gsub('([[:space:]])([A-Z])',' \\L\\2', temp, perl=T))
}

unStuckGeneAnnotationsNames <- function(SpeciesNameInput) {
  temp = sapply(strsplit(SpeciesNameInput, "_"), '[', c(1,2))
  temp = t(temp)
  temp = as.character(temp[,1])
  temp = gsub('([a-z])(?=[A-Z])','\\1 ', temp, perl=T)
  return(gsub('([[:space:]])([A-Z])',' \\L\\2', temp, perl=T))
}

################################################################################
#  ENRICHMENT TOOL                                                             #
#  Extract top CGs/Genes by top Genes OR top CGs                               #
################################################################################
getEnrichmentInputs <- function(cglist = NA, geneMap = NA, genomeVersion = "hg19",
                                topCG = NA, topGenes = NA,
                                rankbywhat = "Z", getGREAT=T,
                                filterSNPs = T, setThreshold = NA) {
  
  myorder = TRUE
  if(rankbywhat == "pvalue") myorder = F
  
  if(filterSNPs == T) {
    cglist = cglist[!grepl("rs", names(cglist))]
  }
  if(is.na(geneMap)[1]) {
    geneMap = readRDS("/Users/caesar/Dropbox (Personal)/MammalianArrayNormalizationTools/geneAnnotations/AnnotationAminHaghani/Latest versions/Human.Homo_sapiens.hg38.Amin.V3.RDS")
  }
  
  if(getGREAT == T & is.numeric(geneMap$geneChr) & grepl("hg", genomeVersion)) {
    geneMap$geneChr = paste0("chr", geneMap$geneChr)
    geneMap$geneChr[geneMap$geneChr=='chr23']='chrX'
    geneMap$geneChr[geneMap$geneChr=='chr24']='chrY'
  } else if(getGREAT == T & is.numeric(geneMap$geneChr) & grepl("mm", genomeVersion)) {
    geneMap$geneChr = paste0("chr", geneMap$geneChr)
    geneMap$geneChr[geneMap$geneChr=='chr20']='chrX'
    geneMap$geneChr[geneMap$geneChr=='chr21']='chrY'
  }
  
  
  input = data.frame(CGid = names(cglist))
  input$cor = cglist; input$CGid = as.character(input$CGid)
  #input$Z = sqrt(n-2)*input$correlation/sqrt(1-input$correlation^2)
  #input$pval = readRDS("allEWAS_pvalue.RDS")[, GetColumn]
  
  #
  input=merge(by='CGid',input,geneMap,all.x=T)
  input=subset(input,!is.na(geneChr))
  #input = input.all[, c("CGid", "cor", "geneChr", "CGstart", "CGend")]
  # order
  input=input[order(abs(input$cor), decreasing = myorder), ]
  # get pos neg background
  pos = input[input$cor > 0, ]
  neg = input[input$cor < 0, ]
  background = input
  
 if(!is.na(topCG) & is.na(topGenes) & !is.na(setThreshold)) {
    
    pos = pos[1:min(topCG, nrow(pos)), ]
    neg = neg[1:min(topCG, nrow(neg)), ]
    
    pvals = pnorm(- abs(pos$cor)) * 2
    pos = pos[pvals <= setThreshold, ]
    pvals = pnorm(- abs(neg$cor)) * 2
    neg = neg[pvals <= setThreshold, ]
    
    print(paste0("Hyper selected: ", nrow(pos), ". Hypo selected: ", nrow(neg)))
    
  } else if(!is.na(topCG) & is.na(topGenes)) {
    
    pos = pos[1:topCG, ]
    neg = neg[1:topCG, ]
      
  } else if(is.na(topCG) & !is.na(topGenes)) {
    for(i in 1:nrow(pos)) {
      geneThreshold = length(unique(pos$ENSEMBL[1:i]))
      if(geneThreshold >= topGenes) {
        geneThreshold = i
        break
      }
    }
    pos = pos[1:geneThreshold, ]
    
    for(i in 1:nrow(neg)) {
      geneThreshold = length(unique(pos$ENSEMBL[1:i]))
      if(geneThreshold >= topGenes) {
        geneThreshold = i
        break
      }
    }
    neg = neg[1:geneThreshold, ]
    
  }
  
  if(getGREAT == T) {
    pos = subset(pos, select=c(geneChr,CGstart,CGend,CGid))
    neg = subset(neg, select=c(geneChr,CGstart,CGend,CGid))
    background = subset(background, select=c(geneChr,CGstart,CGend,CGid))
  }
  
  toreturn = list(pos, neg, background)
  names(toreturn) = c("pos", "neg", "background")
  return(toreturn)
  
}

################################################################################
#  ENRICHMENT TOOL                                                             #
#  One Step GREAT Analysis                                                     #
################################################################################

oneStepGREAT <- function(input = NA, geneMap = NA, setThreshold = NA, topCGNumber = 500,
                         genomeVersion = NA, version = 4) {
  library(rGREAT)
  
  if(!is.list(input)) {
    input = getEnrichmentInputs(cglist = input, geneMap = geneMap, genomeVersion = genomeVersion,
                                topCG = topCGNumber, rankbywhat = "Z",
                                setThreshold = setThreshold)
  }
  extend=c(50,1000)
  
  if(is.na(genomeVersion)) {
    if(version == 4) {
      genomeVersion = "hg38"
    } else if (version == 3) {
      genomeVersion = "hg19"
    } else {
      genomeVersion = "hg19"
    }
  }
  print(paste0("Top positives: ", nrow(input[[1]]), " CGs."))
  
  # & nrow(input[[1]]) >= 5
  if(sum(is.na(input[[1]][1, ])) == 0 & nrow(input[[1]]) >= 3) {
    job = submitGreatJob(input[[1]], bg = input[[3]],
                         species               = genomeVersion,
                         includeCuratedRegDoms = TRUE,
                         rule                  = c("basalPlusExt"),
                         adv_upstream          = 5.0,
                         adv_downstream        = 1.0,
                         adv_span              = extend[1],
                         request_interval = 0,
                         version=as.character(version),
                         max_tries = 10)
    
    ontology.all=availableOntologies(job)
    #availableOntologies(job.pos, category = "Pathway Data")
    #availableOntologies(job.pos, category = "GO")
    output.all1={}
    if(nrow(input[[3]]) <= 1e5) {
      for(k in 1:length(ontology.all)){
        #print(ontology.all[k])
        out0.list = tryCatch(getEnrichmentTables(job,ontology=ontology.all[k], download_by = "tsv"),error=function(e){NULL})
        # out0.list = tryCatch(getEnrichmentTables(job,ontology=ontology.all[k]),error=function(e){NULL})
        if(!is.null(out0.list)){
          if(length(out0.list[[1]]$Ontology) > 0) {
            db0.list=as.list(names(out0.list))
            output<-Map(cbind,Database=db0.list,out0.list)
            output<-do.call('rbind',output)
            output.all1=rbind(output.all1,output)
          }
        }
      }
    } else {
      ## For 300k array
      for(k in 1:length(ontology.all)){
        # print(ontology.all[k])
        #out0.list = tryCatch(getEnrichmentTables(job,ontology=ontology.all[k], download_by = "tsv"),error=function(e){NULL})
        out0.list = tryCatch(getEnrichmentTables(job,ontology=ontology.all[k]),error=function(e){NULL})
        if(!is.null(out0.list)){
          #db0.list=as.list(names(out0.list))
          #output<-Map(cbind,Database=db0.list,out0.list)
          #output<-do.call('rbind',output)
          #output=subset(out0.list[[1]],select=-c(BgRegionNames,BgRegionNames ))
          output=data.frame(Database=ontology.all[k],out0.list[[1]])
          output.all1=rbind(output.all1,output)
        }
      }
      colnames(output.all1)[which(colnames(output.all1) == "Hyper_Raw_PValue")] = "HyperP"
    }
    output.all1$class =  "pos"
  } else {
    output.all1 = NULL
  }
  
  if(sum(is.na(input[[2]][1, ])) == 0 & nrow(input[[2]]) >= 3) {
    job = submitGreatJob(input[[2]], bg = input[[3]],
                         species               = genomeVersion,
                         includeCuratedRegDoms = TRUE,
                         rule                  = c("basalPlusExt"),
                         adv_upstream          = 5.0,
                         adv_downstream        = 1.0,
                         adv_span              = extend[1],
                         request_interval = 0,
                         version= as.character(version),
                         max_tries = 10)
    
    ontology.all=availableOntologies(job)
    #print(ontology.all)
    #availableOntologies(job.pos, category = "Pathway Data")
    #availableOntologies(job.pos, category = "GO")
    output.all={}
    if(nrow(input[[3]]) <= 1e5) {
      for(k in 1:length(ontology.all)){
        out0.list = tryCatch(getEnrichmentTables(job,ontology=ontology.all[k], download_by = "tsv"),error=function(e){NULL})
        if(!is.null(out0.list)){
          if(length(out0.list[[1]]$Ontology) > 0) {
            db0.list=as.list(names(out0.list))
            output<-Map(cbind,Database=db0.list,out0.list)
            output<-do.call('rbind',output)
            output.all=rbind(output.all,output)
          }
        }
        
      }
      
    } else {
      ## For 300k array
      for(k in 1:length(ontology.all)){
        # print(ontology.all[k])
        #out0.list = tryCatch(getEnrichmentTables(job,ontology=ontology.all[k], download_by = "tsv"),error=function(e){NULL})
        out0.list = tryCatch(getEnrichmentTables(job,ontology=ontology.all[k]),error=function(e){NULL})
        if(!is.null(out0.list)){
          #db0.list=as.list(names(out0.list))
          #output<-Map(cbind,Database=db0.list,out0.list)
          #output<-do.call('rbind',output)
          #output=subset(out0.list[[1]],select=-c(BgRegionNames,BgRegionNames ))
          output=data.frame(Database=ontology.all[k],out0.list[[1]])
          output.all=rbind(output.all,output)
        }
      }
      colnames(output.all)[which(colnames(output.all) == "Hyper_Raw_PValue")] = "HyperP"
    }
    output.all$class =  "neg"
    if(!is.null(output.all1)) output.all = rbind(output.all1, output.all)
    
  } else {
    if(is.null(output.all1)) {
      output.all = NULL
    } else {
      output.all = output.all1
    }
  }
  
  # output.all = rbind(output.all1, output.all)
  if(!is.null(output.all)) {
    if(nrow(output.all) > 2) {
      output.all = output.all[order(output.all$HyperP, decreasing = FALSE), ]
      output.all = output.all[output.all$HyperP <= 0.05, ]
    }
  }
  return(output.all)
}

################################################################################
#  Stouffer Summary                                                            #
#  Summarize over tissues                                                      #
################################################################################
getStouffer <- function(mydf, myweights = NA) {
  if(is.na(myweights[1])) {
    myweights = sapply(strsplit(colnames(mydf), "\\.N"), "[", c(2))
    myweights = sqrt(as.numeric(myweights))
  } 
  notna = !is.na(mydf[1, ])
  myweights = myweights[notna]
  if(!length(myweights) == ncol(mydf[, notna])) 
    stop("Number of weights do not match number of columns")
  displayNames = sapply(strsplit(colnames(mydf[, notna]), "Tissue"), '[', c(2))
  displayNames = sapply(strsplit(displayNames, "\\.N"), '[', c(1))
  print(paste0(paste(displayNames, collapse = " "), " columns are being summarized."))
  return(apply(mydf[, notna], 1, function(x) return(sum(x * myweights) / sqrt(sum(myweights^2)))))
}


getGeneDiscription <- function() {
  library(biomaRt)
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  IDs <- c("BRCA2","BRAF")
  genedesc <- getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = IDs, mart =ensembl)
  head(genedesc)
}


################################################################################
#  EWASGWAS Enrichment                                                         #
#                                                                              #
################################################################################
F_enrich<- function(ORDER,aar0,other0,other.name0,cutoff,n.topcpg0,background.hg=background.hg19,
                    mythreshold = NA){
  aar0=subset(aar0,!is.na(aar0$Meta.Z))
  other0=subset(other0,!is.na(other0$p.other))
  other0.gene=other0
  other0=merge(by='HGNC_gene',other0,subset(background.hg,select=c(HGNC_gene,CHR,bp,CGid,HGNC_gene.amin)))
  other0$genomic_region=paste(as.character(other0$HGNC_gene),as.character(other0$CGid),sep='-')
  #length(unique(other0$GeneID))
  #
  background=other0#no chr23 and chr24
  #very important
  ntot=nrow(background)
  #
  aar0=aar0[order(-abs(aar0$Meta.Z)),]
  aar0$genomic_region=paste(as.character(aar0$HGNC_gene),as.character(aar0$CGid),sep='-')
  aar.top=vector(length=2,mode='list')
  aar.top[[1]]=subset(aar0,Meta.Z>0)
  aar.top[[2]]=subset(aar0,Meta.Z<0)
  aar.top[[1]]=aar.top[[1]][1:n.topcpg0,]
  aar.top[[2]]=aar.top[[2]][1:n.topcpg0,]
  if(!is.na(mythreshold)) {
    aar.top[[1]]=aar.top[[1]][aar.top[[1]]$P.value <= mythreshold, ]
    aar.top[[2]]=aar.top[[2]][aar.top[[2]]$P.value <= mythreshold,]
  }
  
  names(aar.top)=c('pos','neg')
  #
  other0=other0[order(other0$p.other),]
  other0.gene=other0.gene[order(other0.gene$p.other),]
  #
  output={}
  for(kk in 1:length(cutoff)){
    n2=round(cutoff[kk]*dim(other0.gene)[1])#top cutoff genes
    other.gene=other0.gene[1:n2,]
    #
    other=subset(other0,HGNC_gene%in%other.gene$HGNC_gene)
    for(t in 1:2){
      index=is.element(aar.top[[t]]$genomic_region,other$genomic_region)
      x=sum(as.numeric(index))
      m=dim(other)[1]# number white balls/other gwas
      n=ntot-m
      k=dim(aar.top[[t]])[1] # numer of draw
      p.enrich=1-sum(dhyper(0:x-1,m=m,n=n,k=k)) #obser x or > x
      if(p.enrich==0){p.enrich=5.0e-17}
      #
      OverlapGenes=paste(unique(aar.top[[t]]$HGNC_gene[index]),collapse =';')
      OverlapGenes.CpG=paste(aar.top[[t]]$CGid[index],collapse =';')
      output0=data.frame(Index=ORDER,class=names(aar.top)[t],cutoff=cutoff[kk],GWAS=other.name0,Overlap=paste(x,m,sep="/"),
                         P=p.enrich, OverlapGenes=OverlapGenes,OverlapGenes.CpG=OverlapGenes.CpG)
      output=rbind(output,output0)
    }
    
  }
  return(output)
}
oneStepGWAS <- function(cglist, topNumber = 500, geneMap, 
                        folderpath = "~/Dropbox (Personal)/HorvathLabCoreMembers/Ake/EnrichmentAnalysis/EWASGWAS/",
                        n.topcpg = 1000, cutoff1 = 0.025, mythreshold = NA) {
  library(dplyr)
  
  background.hg18 = read.csv(paste0(folderpath, "input/Magenta_mammlian_background_hg18.csv"), 
                             stringsAsFactors = F)
  background.hg18 = background.hg18[background.hg18$CGid %in% names(cglist), ]
  background.hg19 = read.csv(paste0(folderpath, "input/Magenta_mammlian_background_hg19.csv"), 
                             stringsAsFactors = F)
  background.hg19 = background.hg19[background.hg19$CGid %in% names(cglist), ]
  
  
  aars.hg18 = read.csv(paste0(folderpath, "input/MAGENTA_All_genomic_region_hg18.csv"), 
                       stringsAsFactors = F)
  aars.hg18 = aars.hg18[aars.hg18$CGid %in% names(cglist), ]
  aars.hg19 = read.csv(paste0(folderpath, "input/MAGENTA_All_genomic_region_hg19.csv"), 
                       stringsAsFactors = F)
  aars.hg19 = aars.hg19[aars.hg19$CGid %in% names(cglist), ]
  
  ## aars = list of all cgs and their gene annotations
  
  aars = data.frame(CGid = names(cglist), Meta.Z = cglist)
  aars$P.value = 2 * pnorm(-abs(aars$Meta.Z)) 
  
  aars.hg18 = aars.hg18[, -which(colnames(aars.hg18) %in% c("Meta.Z", "P.value"))] %>%
    left_join(aars, by = "CGid")
  aars.hg19 = aars.hg19[, -which(colnames(aars.hg19) %in% c("Meta.Z", "P.value"))] %>%
    left_join(aars, by = "CGid")
  
  ## other = gwas annotation, which needs to be fixed from Ake's code
  gwas.anno = read.csv(paste0(folderpath, "GWASDataAnnotation.csv"), stringsAsFactors = F)
  gwas.anno$data = sapply(strsplit(gwas.anno$data, "/"), function(x) return(x[length(x)]))
  output.all={}
  for(k in 1:nrow(gwas.anno)){
    other=read.csv(gzfile(paste0(folderpath, "GWASdata/", gwas.anno$data[k])))
    other.name=gwas.anno$trait.short[k]
    if(gwas.anno$hg[k]=='hg19'){
      #(ORDER           ,aar0,other0,other.name0,cutoff,n.topcpg0,background.hg=background.hg19)
      output0 = F_enrich(gwas.anno$Index[k],aars.hg19,other,other.name,cutoff1,n.topcpg,background.hg=background.hg19,
                         mythreshold = mythreshold)
    } else {
      output0 = F_enrich(gwas.anno$Index[k],aars.hg18,other,other.name,cutoff1,n.topcpg,background.hg=background.hg18,
                         mythreshold = mythreshold)
    }
    output.all=rbind(output.all,output0)
  }
  output.all=merge(by='Index',gwas.anno,output.all)
  output.all=subset(output.all,select=-c(data,hg,trait.short))
  output.all=output.all[order(output.all$Index),]
  #
  
  if(grepl("pos|neg", output.all$class[1])) {
    output.all$P_FDR = NA
    output.all[grepl("pos", output.all$class), "P_FDR"] = p.adjust(output.all[grepl("pos", output.all$class), "P"], method = "fdr")
    output.all[grepl("neg", output.all$class), "P_FDR"] = p.adjust(output.all[grepl("neg", output.all$class), "P"], method = "fdr")
  }
  
  return(output.all)
}
