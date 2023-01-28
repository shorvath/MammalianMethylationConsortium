rm(list=ls())
options(stringsAsFactors = F)
#setwd('')
meth.file=c('simulated_450k.rds')
impute.qc=read.csv('MeasuresTEST_EPICtoMM_v2_450k.csv')
######## load the functions
source("fns_ArrayConverterAlgorithm.R")

######## load the trained model
fitall = readRDS('fit_EPICtoMM_v2_EPIC450k.rds')

########## a quick look at the trained models
length(fitall)
head(names(fitall))
head(fitall[[1]])

###### load your Illumina array data
dat.meth=readRDS(meth.file)
rownames(dat.meth)=dat.meth$ID
dat.meth$ID<-NULL
###### format the array with rownames as the CG names and columns are the samples.

#
cpg.need=unlist(lapply(fitall,'[[',1))
cpg.need=unique(cpg.need)
cpg.need=cpg.need[cpg.need!="(Intercept)"]
#
cpg.miss=is.element(cpg.need,rownames(dat.meth))
table(cpg.miss)
cpg.miss=cpg.need[!cpg.miss]
#
#impute missing cpg by beta=0.5
#
if (length(cpg.miss)>0){
  add=matrix(0.5,nrow=length(cpg.miss),ncol=ncol(dat.meth))
rownames(add)=cpg.miss
colnames(add)=colnames(dat.meth)
dat.meth=rbind(dat.meth,add)
}
#
#re-check any missing cpgs
cpg.miss=is.element(cpg.need,names(dat.meth))
table(cpg.miss)
#
#
########## imputed mammalian CpGs
glment.time<-proc.time()
cglist = names(fitall)
out1 = ArrayConverter_EPICtoMM(cglist,
                               InputEPIC = dat.meth)
dim(out)
saveRDS(out,'myMammalianImputedArray.Rds')
proc.time()-glment.time
#
#QC of imputed Mammalian CpGs based on our test dataset can be found in MeasuresTEST_EPICtoMM_v2_450k.csv
#
#
head(impute.qc)#We suggest that the QC of the imputed CpGs is guided by bicor_test>0.6 & cor_test>0.6