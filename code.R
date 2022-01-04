rm(list=ls())
getwd()
library(readxl)
library(glmnet)
library(openxlsx)
library(caret)
library(survAUC)
library(survival)
library(survminer)
library(mice)
library(plyr)
library(boot)
library(rms)
library(MASS)
# Evaluate the intra-observer repeatability for each radiomic feature
library(irr)
READ11<-read_xlsx(".......xlSX")
READ12<-read_xlsx(".........xlSX")
ICC<-function(x){
  read11<-READ11[ ,x]
  read12<-READ12[ ,x]
  DATA<-cbind(read11,read12)
  ICC<-icc(DATA, model="twoway",type="agreement", unit = "single")
  ICCVALUE<-round(ICC$value,3)
  CI<-paste(round(ICC$lbound,3),sep = "-",round(ICC$ubound,3))
  P<-round(ICC$p.value,3)
  iccvalues<-data.frame( 
    'ICCVALUE'=ICCVALUE,
    '95 CL'=CI,
    'Pvalue'=P)
  return(iccvalues)
}
ICC()

library(plyr)

iccdata<-lapply(c(2:2+n),ICC)
iccdata<-ldply(iccdata,data.frame)

iccdata$feature<-NA

iccdata$feature<-names(READ11[-c(1)])


iccdata_good1<-iccdata[iccdata$ICCVALUE>0.75, ]
# The intra-observer repeatability for each radiomic feature was evaluated in the similar method
iccdata_good2<-iccdata2[iccdata2$ICCVALUE>0.75, ]
# Repeatable features(Features with intraclass correlation coefficient >0.75 in both intra- and inter-observer repeatability tests)
good_feature<-intersect(iccdata_good1$feature,iccdata_good2$feature)

#  feature harmonization

library(neuroCombat)
library(matrixStats)


#Batch variable for the scanner 
batch<-as.numeric(TAL_data$mechine_classification)
v = as.vector(unlist(batch))
DATA<-as.matrix(radiomics_data,rownames.force = NA)
colnames(DATA)<-NULL

data.harmonized <- neuroCombat(dat=DATA, batch=v,eb=FALSE,mean.only=TRUE )

data.harmonized$dat.combat

write.xlsx(data.harmonized$dat.combat,"data.harmonized.xlsx") 


#normalization of  features
set.seed()
sampling_vector<-createDataPartition(TAL_data$follow_outcome,p=0.7,list=F)
sampling_vector<-as.vector(sampling_vector)
traindata<-TAL_data[sampling_vector, ]
trainlable<-TAL_data$follow_outcome[sampling_vector]
testdata<-TAL_data[-sampling_vector, ]
testlable<-TAL_data$follow_outcome[-sampling_vector]

train_numerber<-traindata[ ,c(x:x+n)]
test_numerber<-testdata[ ,c(x:x+n)]
pp_unit<-preProcess(train_numerber,method = c("center","scale"))

train_norm<-predict(pp_unit,train_numerber)
test_norm<-predict(pp_unit,test_numerber)
traindata[ ,c(x:x+n)]<-train_norm
testdata[ ,c(x:x+n)]<-test_norm
# The correlation coefficient between each pair of features

library(caret)
library(matrixStats)

# calculate correlation matrix
correlationMatrix <- cor()
# summarize the correlation matrix
print(correlationMatrix)
# find attributes that are highly corrected 
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.9,names = TRUE)


# The least absolute shrinkage and selection operator (LASSO) penalized Cox-model 
train_mat<-as.matrix(traindata[ ,c(x:x+m)])
test_mat<-as.matrix(testdata[ ,c(x:x+m)])
basesurv<-Surv(time = traindata$OS_m,event=traindata$follow_outcome)

lambdas<-10^seq(-0.1,-4,length=250)
models_lasso<-glmnet(train_mat,basesurv, alpha=1,lambda = lambdas,family = "cox",standardize=F)

plot(models_lasso,xvar="lambda",main="Lasso\n")

lasso_cv<-cv.glmnet(train_mat,basesurv, family = "cox",alpha=1,lambda = lambdas)
lambda_lasso<-lasso_cv$lambda.min
lambda_lasso
plot(lasso_cv)
g.coef <- as.matrix(coef(lasso_cv, s = lambda_lasso))   
c<-cbind(g.coef,colnames(traindata[ ,c(x:x+m)]))
xishu<-c[(which(c[ ,1]!= 0)), ] 

predict_train <- predict(models_lasso,s=lambda_lasso,type = "response",newx=train_mat)
predict_test<-predict(models_lasso,s=lambda_lasso,type = "response",newx=test_mat)


traindata$radscore<-NA
traindata$radscore<-predict_train[ ,1] 

testdata$radscore<-NA
testdata$radscore<-predict_test[ ,1] 

# The radscore was converted into dichotomized radiomic signature


res.cut <- surv_cutpoint(traindata,time ="OS_m",event="follow_outcome",
                         variables = "radscore")

plot(res.cut, "radscore", palette = "npg")

traindata<-within(traindata,{
  radscore_risk<-NA
  radscore_risk<-as.factor(ifelse((radscore>=1.01),1,0))
})
testdata<-within(testdata,{
  radscore_risk<-NA
  radscore_risk<-as.factor(ifelse((radscore>=1.01),1,0))
})


# The clinicopathologic-radiographic (CPR) model

Unicox_c<-function(x){
  FML<-as.formula(paste0('basesurv~',x))
  Ecox<-coxph(FML,data=traindata)
  Esum<-summary(Ecox)
  HR<-round(Esum$coefficients[ ,2],3)
  P<-round(Esum$coefficients[ ,5],3)
  CI<-paste0(round(Esum$conf.int[ ,3:4],3),collapse = "-")
  unicox<-data.frame('charateristics'=x,
                     'Hazard Ratio'= HR,
                     'CI95'= CI,
                     'P.value'= P)
  return(unicox)
}

Varnames<-colnames(traindata)[....]
unicox_list_c<-lapply(Varnames,Unicox_c)
unicox_list_c<-ldply(unicox_list_c,data.frame)


library(MASS)


FML<-as.formula(paste0('basesurv~',paste0(unicox_list_c[unicox_list_c$P.value<0.1, ]$charateristics,collapse = "+")))

mod.none<-Multicox<-coxph(basesurv~1,data=traindata)
mod.full<-Multicox<-coxph(FML,data=traindata)

lm_back<-stepAIC(mod.full,direction = "backward",k=3.8415)
Multicox_cpr<-coxph(basesurv~aa+bb+cc+dd,data=traindata,x=T,y=T)

#The proportional hazards assumption for each variable
Unicoxph<-function(x){
  FML<-as.formula(paste0('basesurv~',x))
  Ecox<-coxph(FML,data=traindata)
  Esum<-summary(Ecox)
  ph <- cox.zph(Ecox)
  return(ph)
}
Varnames<-colnames(traindata)[....]
phph<-lapply(Varnames,Unicoxph)
phph


#  A combined model integrating the radiomic signature with CPR factors
Multicox_cpr_r<-coxph(basesurv~radscore_risk+aa+bb+cc+dd,data=traindata,x=T,y=T)


