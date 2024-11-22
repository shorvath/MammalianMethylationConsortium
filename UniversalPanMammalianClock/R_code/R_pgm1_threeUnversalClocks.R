#
#Universal mammalian clock2 of age
#
#author for developing clock2: Ake T. Lu
#author for developing clock1 and clock3: Zhe Fei
#
rm(list=ls())
options(stringsAsFactors = F)
#
myinput.list=readRDS('mydata_GitHub.Rds')
# The following labels the coefficient values of the three universal clocks
beta.name=c('beta_clock1','beta_clock2','beta_clock3')
y.name=c('Y.pred1','Y.pred2','Y.pred3')
age.name=c('DNAmAgePanMammalianClock1','DNAmAgePanMammalianClock2','DNAmAgePanMammalianClock3')

#clock2
F2_antitrans_clock2<-function(y,y.maxAge,y.gestation,const=1){
  x0=const*exp(-exp(-1*y))
  x1=x0*(y.maxAge+y.gestation)
  x=x1-y.gestation
  x
}
#
#clock3
#
F1_logli <- function(age1, m1, m2 = m1, c1=1){
  ifelse(age1 >= m1, (age1-m1)/m2 , c1*log((age1-m1)/m2/c1 +1) )
}
#RelativeAdultAge
F2_revtrsf_clock3 <- function(y.pred, m1, m2 = m1, c1=1){
  ifelse(y.pred<0, (exp(y.pred/c1)-1)*m2*c1 + m1, y.pred*m2+m1 )
}

############# 
# The `loglifn` function shows how to calculate m1 for the transformation
# It is the `a_Logli` in the function

F3_loglifn = function(dat1,b1=1,max_tage = 4,
                      c1=5, c2 = 0.38, c0=0){
  n=nrow(dat1)
  
  age1 = (dat1$maxAge+dat1$GestationTimeInYears)/(dat1$averagedMaturity.yrs+dat1$GestationTimeInYears)
  
  a1 = age1/(1+max_tage)
  dat1$a1_Logli = a1 #x/m1 in manuscript
  
  a2 = (dat1$GestationTimeInYears + c0)/(dat1$averagedMaturity.yrs) 
  dat1$a_Logli = a_Logli = c1*a2^c2
  #m=5*(G/ASM)^0.38 from regression analysis/formula(7)
  
  
  x = dat1$Age + dat1$GestationTimeInYears
  t2 = dat1$averagedMaturity.yrs*b1 + dat1$GestationTimeInYears
  x2 = x/t2 #### log(x/t2)
  y = F1_logli(x2, a_Logli, a_Logli)
  
  dat1$LogliAge <- y
  return(dat1)
}
#
#(1)generate variable HighmaxAge
#
names(myinput.list)
info=myinput.list[[1]]#The only required variables are SpeciesLatinName and Basename 
anage=myinput.list[[3]]
info=merge(by='SpeciesLatinName',info,subset(anage,select=c(SpeciesLatinName,GestationTimeInYears,
                                                             averagedMaturity.yrs,maxAge)))
head(info)
#Description for mymax=1.3
#We were concerned that the uneven evidence surrounding the maximum age of different species 
#could bias our analysis. While billions of people have been evaluated for estimating 
#the maximum age of humans (122.5 years) or mice (4 years), 
#the same cannot be said for any other species. 
#To address this concern, we made the following assumption: 
#the true maximum age is 30% higher than that reported in AnAge
#for all species except for humans and mice (Mus musculus). 
#Therefore, we multiplied the reported maximum lifespan of non-human or non-mouse species by 1.3. 
#Our predictive models turn out to be highly robust with respect to this assumption.
MYMAX=1.3
info$HighmaxAge=MYMAX*info$maxAge
info$HighmaxAge[info$SpeciesLatinName=='Homo sapiens']=info$maxAge[info$SpeciesLatinName=='Homo sapiens']
info$HighmaxAge[info$SpeciesLatinName=='Mus musculus']=info$maxAge[info$SpeciesLatinName=='Mus musculus']
#
#(2)merge info and metharray beta values
#
glmnet.list=myinput.list[[4]]#The three universal clock prediction models
mycpgs=c(glmnet.list[[1]]$var,glmnet.list[[2]]$var,glmnet.list[[3]]$var)
mycpgs=unique(mycpgs)
mycpgs=mycpgs[mycpgs!='Intercept']
#dat.meth0: number of Mammalian array CpGs (n=37554) x [number of samples+1]
#
dat.meth0=myinput.list[[2]]
dat.meth0=subset(dat.meth0,CGid%in%mycpgs)#only keep the CpGs in the three clocks
dat.meth=t(dat.meth0[,-c(1)])
colnames(dat.meth)=dat.meth0$CGid
dat.meth=data.frame(Basename=colnames(dat.meth0)[-c(1)],dat.meth)
dat.meth$Intercept=1
#
info=merge(by='Basename',info,dat.meth)
#
#predict RelativeAge
#

for(k in 1:3){
glmnet=glmnet.list[[k]]
glmnet$beta=glmnet[,beta.name[k]]
#glmnet$var[1]=ifelse(glmnet$var[1]=="(Intercept)",'Intercept',glmnet$var[1])
temp=as.matrix(subset(info,select=as.character(glmnet$var)))
info[,y.name[k]]=as.numeric(as.matrix(subset(info,select=as.character(glmnet$var)))%*%glmnet$beta)
}
#(1) Clock 1
info[,age.name[1]]=exp(info[,y.name[k]])-2
#(2) Clock 2
info$DNAmRelativeAge=exp(-exp(-1*info[,y.name[2]]))
info[,age.name[2]]= F2_antitrans_clock2(info[,y.name[2]],info$HighmaxAge,info$GestationTimeInYears,const=1)
#(3) Clock 3
info=F3_loglifn(info)#to compute m estimate for tuning point in the log-linear transformation
info$m1=info$a_Logli
info$DNAmRelativeAdultAge=F2_revtrsf_clock3(info[,y.name[3]], info$m1)
info[,age.name[3]]<-
  info$DNAmRelativeAdultAge *(info$averagedMaturity.yrs + info$GestationTimeInYears) -info$GestationTimeInYears
#
output=subset(info,select=c('Basename','SpeciesLatinName','MammalNumberHorvath','Age','Tissue','DNAmRelativeAge','DNAmRelativeAdultAge',age.name))
#
x=output$Age
y=output$DNAmAgePanMammalianClock2
cor0=round(cor(x,y),2)
mae0=round(median(abs(x-y)),2)
plot(x,y,xlab='Age',ylab='Clock2: DNAm Age',
     main=paste0('Bottlenose dolphin (n=',length(x),')\n','cor=',cor0,
                 ', MAE=',mae0),
     col='white',xlim=c(0,max(max(x),max(y))+1),
     ylim=c(0,max(max(x),max(y))+1))
text(x,y,label=output$MammalNumberHorvath,col='blue')
#

