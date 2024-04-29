#Dependancies------
library(dplyr) # data wrangling
library(tidyr) # data wrangling
library(readxl) # read in excel files
library(car) # Anova types
library(visreg) # data visualizations
library(glmulti) # model selection
library(lme4) # random mixed effect models
library(lmerTest) # test random effects
library(stargazer) # publishable tables
library(sjPlot) # publishable tables
library(stringr) # paste characters
library(modelsummary) #goodness of fit in LME model
library(AICcmodavg) # AIC model averaging and comparison
library(khroma) # color blind friendly color palettes 
library(MuMIn) # LME mR2 and cR2

sigfig <- function(vec, n=5){ 
  formatC(signif(vec,digits=n), digits=n,format="fg", flag="#")}

lmer.glmulti<-function(formula,data,random="",...) {
  newf <- formula
  newf[[3]] <- substitute(f+r,
                          list(f=newf[[3]],
                               r=reformulate(random)[[2]]))
  lmer(newf,data=data,
       REML=FALSE,...)
}

glmulti_topfnct<-function(res=res){
  top <- weightable(res)
  top <- top[top$aic <= min(top$aic) + 10000,]
  print(top$model[1])
  return(top)
  remove(res,df,response,m1)
}

glmulti_topfnct2<-function(res=res){
  top <- weightable(res)
  top <- top[top$aic <= min(top$aic) + 2,]
  print(top$model[1])
  return(top)
  remove(res,df,response,m1)
}

wdir="~/Dropbox/KU/Dimensions/SoilRespiration"

#Load Data-----

SoilRespFile="RFiles/Manuscript_Submission/SoilRespiration2020-2023.csv"
SR20.23<-read.csv(file=paste(wdir,SoilRespFile,sep="/"))
SR20.23<-SR20.23[,-match("X",names(SR20.23))]
str(SR20.23)
vect<-match(c("MONTH","RAINTRT","DISPERSION","BLOCK","SUBBLOCK","DISPERSION.FOURLEVELS"),names(SR20.23))
for (i in 1:NROW(vect)){SR20.23[,vect[i]]<-as.factor(SR20.23[,vect[i]])}
remove(vect)

SR20.23.Calc<-subset(SR20.23,Moisture.Calc>=0)

SR20.23.Calc$ScaledMoisture<-scale(SR20.23.Calc$Moisture.Calc)
SR20.23.Calc$ScaledMoisture2<-scale(SR20.23.Calc$Moisture.Calc2)
SR20.23.Calc$YEARMONTH<-as.factor(paste(SR20.23.Calc$YEAR,SR20.23.Calc$MONTH, sep=" "))

SR20.23.anti<-subset(SR20.23,VWC.perc>0)
SR20.23.anti$ScaledMoisture<-scale(SR20.23.anti$VWC.perc)
SR20.23.anti$ScaledMoisture2<-scale(SR20.23.anti$VWC.perc2)
SR20.23.anti$YEARMONTH<-as.factor(paste(SR20.23.anti$YEAR,SR20.23.anti$MONTH, sep=" "))

DF5cm<-subset(SR20.23.Calc,YEARMONTH=="2020 10"|YEARMONTH=="2020 5"|YEARMONTH=="2020 6"|YEARMONTH=="2020 7"|YEARMONTH=="2020 8"|YEARMONTH=="2020 9"|YEARMONTH=="2022 5"|YEARMONTH=="2022 6")
DF12cm<-subset(SR20.23.anti,YEARMONTH=="2021 10"|YEARMONTH=="2021 5"|YEARMONTH=="2021 6"|YEARMONTH=="2021 7"|YEARMONTH=="2021 8"|YEARMONTH=="2021 9"|YEARMONTH=="2022 7"|YEARMONTH=="2022 8"|YEARMONTH=="2022 9"|YEARMONTH=="2023 10"|YEARMONTH=="2023 5"|YEARMONTH=="2023 6"|YEARMONTH=="2023 7"|YEARMONTH=="2023 8"|YEARMONTH=="2023 9"|YEARMONTH=="2021 10"|YEARMONTH=="2021 6"|YEARMONTH=="2021 7"|YEARMONTH=="2021 8"|YEARMONTH=="2021 9")

names(DF5cm)==names(DF12cm)

SR20.23.scale<-rbind(DF5cm,DF12cm)

SR20.23.9<-subset(SR20.23,MONTH==9)
SR20.22.9<-subset(SR20.23.9,YEAR!=2021)
SR20.22.9<-subset(SR20.22.9,YEAR!=2023)

SR20.23.Calc<-subset(SR20.23,Moisture.Calc>=0)

#EDA-----
#Soil Respiration
hist(SR20.23$FLUX.mgCO2C.m2.h)
hist(log10(SR20.23$FLUX.mgCO2C.m2.h))

#Soil Moisture and Temperature
hist(SR20.23$Moisture.Calc)
hist(SR20.23$VWC.perc)
hist(SR20.23$Temp.C)

#Microbial and Root Biomass
hist(SR20.23$DryWt.g.m2)
hist(log(SR20.23$DryWt.g.m2))
hist(SR20.23$MB.mg.gdw)
hist(log(SR20.23$MB.mg.gdw))
       
#Covariate Models-----------
##AICmodelavg----
covar.ml<-lmer(log10(FLUX.mgCO2C.m2.h) ~ scale(VWC.perc) + scale(VWC.perc2) + Temp.C 
            + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
            , data=SR20.23,REML=F)

covar.1.ml<-lmer(log10(FLUX.mgCO2C.m2.h) ~ scale(VWC.perc) + scale(VWC.perc2) + Temp.C 
              + MONTH + YEAR
              + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK)
              , data=SR20.23,REML=F)

covar.2.ml<-lmer(log10(FLUX.mgCO2C.m2.h) ~ scale(VWC.perc) + scale(VWC.perc2) 
              + MONTH + YEAR
              + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK)
              , data=SR20.23,REML=F)

covar.3.ml<-lmer(log10(FLUX.mgCO2C.m2.h) ~ 
                Temp.C 
              + MONTH + YEAR
              + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK)
              , data=SR20.23,REML=F)

covar.4.ml<-lmer(log10(FLUX.mgCO2C.m2.h) ~ 
                Temp.C 
              + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK)
              , data=SR20.23,REML=F)

covar.5.ml<-lmer(log10(FLUX.mgCO2C.m2.h) ~ scale(VWC.perc) + scale(VWC.perc2) + 
                + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK)
              , data=SR20.23,REML=F)

AIC(covar.ml,covar.1.ml,covar.2.ml,covar.3.ml,covar.4.ml,covar.5.ml)

CovariateModels12<-aictab(Cand.models <- list(
  "12cm Moisture + 12cm Moisture2 + Temp" = covar.ml,
  "12cm Moisture + 12cm Moisture2 + Temp + Month + Year" = covar.1.ml,
  "12cm Moisture + 12cm Moisture2 + Month + Year" = covar.2.ml,
  "Temp + Month + Year" = covar.3.ml,
  "Temp" = covar.4.ml,
  "12cm Moisture + 12cm Moisture2" = covar.5.ml
))

CovariateModels12

Moisture5cm.ml<-lmer(log10(FLUX.mgCO2C.m2.h)~
                       Moisture.Calc + Moisture.Calc2 +
                       + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
                     , data=SR20.23.Calc,REML=F)

Moisture5cm.1.ml<-lmer(log10(FLUX.mgCO2C.m2.h)~
                         Moisture.Calc +  Moisture.Calc2 + MONTH + YEAR 
                       + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
                       , data=SR20.23.Calc,REML=F)

Moisture5cm.2.ml<-lmer(log10(FLUX.mgCO2C.m2.h)~
                         Moisture.Calc + Moisture.Calc2 + Temp.C +
                         + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
                       , data=SR20.23.Calc,REML=F)

Moisture5cm.3.ml<-lmer(log10(FLUX.mgCO2C.m2.h)~
                         Moisture.Calc + Moisture.Calc2 + Temp.C + MONTH + YEAR 
                       + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
                       , data=SR20.23.Calc,REML=F)

Moisture5cm.4.ml<-lmer(log10(FLUX.mgCO2C.m2.h)~
                         Temp.C + MONTH + YEAR 
                       + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
                       , data=SR20.23.Calc,REML=F)

Moisture5cm.5.ml<-lmer(log10(FLUX.mgCO2C.m2.h)~
                         Temp.C
                       + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
                       , data=SR20.23.Calc,REML=F)

CovariateModels05<-aictab(Cand.models <- list(
  "5cm Moisture + 5cm Moisture2 + Temp" = Moisture5cm.2.ml,
  "5cm Moisture + 5cm Moisture2 + Temp + Month + Year" = Moisture5cm.3.ml,
  "5cm Moisture + 5cm Moisture2 + Month + Year" = Moisture5cm.1.ml,
  "Temp + Month + Year" = Moisture5cm.4.ml,
  "Temp" = Moisture5cm.5.ml,
  "5cm Moisture + 5cm Moisture2" = Moisture5cm.ml
))

CovariateModels05

##REML models----

anova(covar<-lmer(log10(FLUX.mgCO2C.m2.h) ~ scale(VWC.perc) + scale(VWC.perc2) + Temp.C 
               + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
               , data=SR20.23,REML=T))

anova(covar.1<-lmer(log10(FLUX.mgCO2C.m2.h) ~ scale(VWC.perc) + scale(VWC.perc2) + Temp.C 
                 + MONTH + YEAR
                 + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK)
                 , data=SR20.23,REML=T))

anova(covar.2<-lmer(log10(FLUX.mgCO2C.m2.h) ~ scale(VWC.perc) + scale(VWC.perc2) 
                 + MONTH + YEAR
                 + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK)
                 , data=SR20.23,REML=T))

anova(covar.3<-lmer(log10(FLUX.mgCO2C.m2.h) ~ 
                   Temp.C 
                 + MONTH + YEAR
                 + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK)
                 , data=SR20.23,REML=T))

anova(covar.4<-lmer(log10(FLUX.mgCO2C.m2.h) ~ 
                   Temp.C 
                 + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK)
                 , data=SR20.23,REML=T))

anova(covar.5<-lmer(log10(FLUX.mgCO2C.m2.h) ~ scale(VWC.perc) + scale(VWC.perc2) + 

                   + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK)
                 , data=SR20.23,REML=T))


ova(Moisture5cm<-lmer(log10(FLUX.mgCO2C.m2.h)~
                             Moisture.Calc + Moisture.Calc2 +
                             + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
                           , data=SR20.23.Calc,REML=T))

anova(Moisture5cm.1<-lmer(log10(FLUX.mgCO2C.m2.h)~
                               Moisture.Calc +  Moisture.Calc2 + MONTH + YEAR 
                             + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
                             , data=SR20.23.Calc,REML=T))

anova(Moisture5cm.2<-lmer(log10(FLUX.mgCO2C.m2.h)~
                               Moisture.Calc + Moisture.Calc2 + Temp.C +
                               + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
                             , data=SR20.23.Calc,REML=T))

anova(Moisture5cm.3<-lmer(log10(FLUX.mgCO2C.m2.h)~
                               Moisture.Calc + Moisture.Calc2 + Temp.C + MONTH + YEAR 
                             + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
                             , data=SR20.23.Calc,REML=T))

anova(Moisture5cm.4<-lmer(log10(FLUX.mgCO2C.m2.h)~
                               Temp.C + MONTH + YEAR 
                             + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
                             , data=SR20.23.Calc,REML=T))

anova(Moisture5cm.5<-lmer(log10(FLUX.mgCO2C.m2.h)~
                               Temp.C
                             + (1|PLOTNO) + (1|RAINTRT:SUBBLOCK) + (1|SUBBLOCK) 
                             , data=SR20.23.Calc,REML=T))


#Diversity Models-----
##SR: SM 12cm-----

model.2<-lmer(log10(FLUX.mgCO2C.m2.h)~
                YEAR * RAINTRT  * DISPERSION.FOURLEVELS * SPPNO * MONTH
              + scale(VWC.perc) + scale(VWC.perc2) + PROPTOTPLOT
              + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) + (1|PLOTNO)
              , data=SR20.23,REML=T)

plot(model.2)
anova(model.2)
rand(model.2)
r.squaredGLMM(model.2)


##SR: SM 5cm------

model.3<-lmer(log10(FLUX.mgCO2C.m2.h)~
                YEAR * RAINTRT * DISPERSION.FOURLEVELS * SPPNO * MONTH
              + scale(Moisture.Calc) + scale(Moisture.Calc2) + PROPTOTPLOT                                      
              + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) + (1|PLOTNO)
              , data=SR20.23.Calc, REML=T)

plot(model.3)
anova(model.3)
rand(model.3)
r.squaredGLMM(model.3)

#SR GlMulti-------

SR20.23$log10SR<-log10(SR20.23$FLUX.mgCO2C.m2.h)

##Running GLMULTI RSOIL-------

df=SR20.23
response=df$log10SR
m1 <- lmer(response ~  RAINTRT  * DISPERSION.FOURLEVELS * SPPNO * YEAR * MONTH
           + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) + (1|PLOTNO)  ,data=df)
res <- glmulti(formula(m1,fixed.only=TRUE),
               random="+(1|PLOTNO)+(1|SUBBLOCK)+(1|RAINTRT:SUBBLOCK)",
               data=df,method="d",
               #deltaM=0.5,
               fitfunc=lmer.glmulti,
               intercept=TRUE,marginality=FALSE,level=2)
res <- glmulti(formula(m1,fixed.only=TRUE),
               random="+(1|PLOTNO)+(1|SUBBLOCK)+(1|RAINTRT:SUBBLOCK)",
               data=df,method="h",
               #deltaM=0.5,
               fitfunc=lmer.glmulti,
               intercept=TRUE,marginality=FALSE,level=2,
               confsetsize = 2916)
SRmodelMNTHYR<-glmulti_topfnct(res=res)
FolderName="CSVoutputfiles"
FileName="SRlmerMNTHYR.csv"
write.csv(SRmodelMNTHYR,file=paste(wdir,FolderName,FileName,sep="/"))


##Loading Model Suites-------

df=SR20.23
response=df$log10SR

FolderName="CSVoutputfiles"
FileName="SRlmerMNTHYR.csv"
SRmodelMNTHYR<-read.csv(file=paste(wdir,FolderName,FileName,sep="/"))
SRmodelMNTHYR<-SRmodelMNTHYR[,-1]
topSRmodelMNTHYR<-SRmodelMNTHYR[1,]
ParaWts.SRMNTHYR<-SRmodelMNTHYR


##Top SR LMER Models-------

m7<-lmer(paste(topSRmodelMNTHYR[1],"+ (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) + (1|PLOTNO)",sep=" "),data=df)
anova(m7)

m7.cvr<-lmer(response ~ RAINTRT + DISPERSION.FOURLEVELS + MONTH + SPPNO + YEAR + MONTH:RAINTRT 
             + MONTH:DISPERSION.FOURLEVELS + RAINTRT:SPPNO + RAINTRT:YEAR 
             + DISPERSION.FOURLEVELS:SPPNO + MONTH:SPPNO + MONTH:YEAR 
             + PROPTOTPLOT
             + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) + (1|PLOTNO),data=df)

m7.cvrmst<-lmer(response ~ RAINTRT + DISPERSION.FOURLEVELS + MONTH + SPPNO + YEAR + MONTH:RAINTRT 
                + MONTH:DISPERSION.FOURLEVELS + RAINTRT:SPPNO + RAINTRT:YEAR 
                + DISPERSION.FOURLEVELS:SPPNO + MONTH:SPPNO + MONTH:YEAR 
                + PROPTOTPLOT + scale(VWC.perc) + scale(VWC.perc2) 
                + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) + (1|PLOTNO),data=df)

m7.cvrmst1<-lmer(response ~ RAINTRT + DISPERSION.FOURLEVELS + MONTH + SPPNO + YEAR + MONTH:RAINTRT 
                 + MONTH:DISPERSION.FOURLEVELS + RAINTRT:SPPNO + RAINTRT:YEAR 
                 + DISPERSION.FOURLEVELS:SPPNO + MONTH:SPPNO + MONTH:YEAR 
                 + PROPTOTPLOT + scale(VWC.perc) 
                 + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) + (1|PLOTNO),data=df)

m7.cvrmst2<-lmer(response ~ RAINTRT + DISPERSION.FOURLEVELS + MONTH + SPPNO + YEAR + MONTH:RAINTRT 
                 + MONTH:DISPERSION.FOURLEVELS + RAINTRT:SPPNO + RAINTRT:YEAR 
                 + DISPERSION.FOURLEVELS:SPPNO + MONTH:SPPNO + MONTH:YEAR 
                 + PROPTOTPLOT 
                 + scale(VWC.perc2) 
                 + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) + (1|PLOTNO),data=df)

m7ml<-lmer(paste(topSRmodelMNTHYR[1],"+ (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) + (1|PLOTNO)",sep=" "),data=df,REML=F)

m7.cvrml<-lmer(response ~ RAINTRT + DISPERSION.FOURLEVELS + MONTH + SPPNO + YEAR + MONTH:RAINTRT 
               + MONTH:DISPERSION.FOURLEVELS + RAINTRT:SPPNO + RAINTRT:YEAR 
               + DISPERSION.FOURLEVELS:SPPNO + MONTH:SPPNO + MONTH:YEAR 
               + PROPTOTPLOT
               + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) + (1|PLOTNO),data=df,REML=F)

m7.cvrmstml<-lmer(response ~ RAINTRT + DISPERSION.FOURLEVELS + MONTH + SPPNO + YEAR + MONTH:RAINTRT 
                  + MONTH:DISPERSION.FOURLEVELS + RAINTRT:SPPNO + RAINTRT:YEAR 
                  + DISPERSION.FOURLEVELS:SPPNO + MONTH:SPPNO + MONTH:YEAR 
                  + PROPTOTPLOT + scale(VWC.perc) + scale(VWC.perc2) 
                  + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) + (1|PLOTNO),data=df,REML=F)

m7.cvrmst1ml<-lmer(response ~ RAINTRT + DISPERSION.FOURLEVELS + MONTH + SPPNO + YEAR + MONTH:RAINTRT 
                   + MONTH:DISPERSION.FOURLEVELS + RAINTRT:SPPNO + RAINTRT:YEAR 
                   + DISPERSION.FOURLEVELS:SPPNO + MONTH:SPPNO + MONTH:YEAR 
                   + PROPTOTPLOT + scale(VWC.perc) 
                   + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) + (1|PLOTNO),data=df,REML=F)

m7.cvrmst2ml<-lmer(response ~ RAINTRT + DISPERSION.FOURLEVELS + MONTH + SPPNO + YEAR + MONTH:RAINTRT 
                   + MONTH:DISPERSION.FOURLEVELS + RAINTRT:SPPNO + RAINTRT:YEAR 
                   + DISPERSION.FOURLEVELS:SPPNO + MONTH:SPPNO + MONTH:YEAR 
                   + PROPTOTPLOT + scale(VWC.perc2) 
                   + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) + (1|PLOTNO),data=df,REML=F)


plot(m7)
plot(m7.cvr)
plot(m7.cvrmst)

anova(m7)
anova(m7.cvr)
anova(m7.cvrmst)

rand(m7)
rand(m7.cvr)
rand(m7.cvrmst)

r.squaredGLMM(m7)
r.squaredGLMM(m7.cvr)
r.squaredGLMM(m7.cvrmst)

aictab(Cand.models <- list(
  "Exp Desgn" = m7ml,
  "Exp Desgn + Cvr " = m7.cvrml,
  "Exp Desgn + Cvr + Mst" = m7.cvrmst1ml,
  "Exp Desgn + Cvr + Mst2" = m7.cvrmst2ml,
  "Exp Desgn + Cvr + Mst + Mst2" = m7.cvrmstml
))


#SM depth Comparison LMER Models----

SR20.23$YEARMONTHPLOT<-paste(SR20.23$YEAR,SR20.23$MONTH,SR20.23$PLOTNO)
SM12<-SR20.23[,match(c("YEARMONTHPLOT","SPPNO","DISPERSION.FOURLEVELS",
                       "RAINTRT","SUBBLOCK","PLOTNO","MONTH","YEAR","VWC.perc"),names(SR20.23))]
SM12$Depth<-as.factor("Depth 12cm")
names(SM12)[9]<-"SM"
SM05<-SR20.23[,match(c("YEARMONTHPLOT","SPPNO","DISPERSION.FOURLEVELS",
                       "RAINTRT","SUBBLOCK","PLOTNO","MONTH","YEAR","Moisture.Calc"),names(SR20.23))]
SM05$Depth<-as.factor("Depth 5cm")
names(SM05)[9]<-"SM"

SMdepth<-rbind(SM12,SM05)
SMdepth<-subset(SMdepth,SM>0)

anova(mod<-lmer(SM~SPPNO*DISPERSION.FOURLEVELS*RAINTRT*Depth
                +(1|PLOTNO)+(1|SUBBLOCK)+(1|RAINTRT:SUBBLOCK)+(1|MONTH)+(1|YEAR)+(1|MONTH:YEAR), 
                data=SMdepth))

plot(mod)
anova(mod)
rand(mod)
r.squaredGLMM(mod)

SMdepth12<-subset(SMdepth,Depth=="Depth 12cm")
SMdepth12<-subset(SMdepth12,YEAR!=2020)

anova(mod2<-lmer(SM~SPPNO
                 *RAINTRT
                 *YEAR
                 +(1|PLOTNO)+(1|SUBBLOCK)+(1|RAINTRT:SUBBLOCK)+(1|MONTH)
                 +(1|DISPERSION.FOURLEVELS) 
                 ,data=SMdepth12))

plot(mod2)
anova(mod2)
rand(mod2)
r.squaredGLMM(mod2)


#Microbial Biomass VS Root Biomass------

SR20<-subset(SR20.23,YEAR==2020)
SR22<-subset(SR20.23,YEAR==2022)

SR20<-subset(SR20,SPPNO!=0)
SR22<-subset(SR22,SPPNO!=0)

SR20.22<-subset(SR20.23,YEAR==2020|YEAR==2022)

SR20.9<-subset(SR20.22.9,YEAR==2020)
SR22.9<-subset(SR20.22.9,YEAR==2022)

SR20.9$ScaledMoisture<-scale(SR20.9$Moisture.Calc)
SR20.9$ScaledMoisture2<-scale(SR20.9$Moisture.Calc2)
SR22.9$ScaledMoisture<-scale(SR22.9$VWC.perc)
SR22.9$ScaledMoisture2<-scale(SR22.9$VWC.perc2)

SR20.22.9<-rbind(SR20.9,SR22.9)

SR.Sept.1reml<-lmer(log10(FLUX.mgCO2C.m2.h)~
                      RAINTRT * DISPERSION.FOURLEVELS * SPPNO
                    
                    + scale(log(MB.mg.gdw))
                    + scale(PROPTOTPLOT)  
                    + scale(lnRootBiomass)
                    + ScaledMoisture                              
                    + ScaledMoisture2
                    
                    + (1|PLOTNO) + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) 
                    , data=SR20.22.9,REML=T)

SR.Sept.2reml<-lmer(log10(FLUX.mgCO2C.m2.h)~
                      RAINTRT * DISPERSION.FOURLEVELS * SPPNO
                    
                    + scale(log(MB.mg.gdw))
                    #+ scale(PROPTOTPLOT)  
                    + scale(lnRootBiomass)
                    + ScaledMoisture                              
                    + ScaledMoisture2
                    
                    + (1|PLOTNO) + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) 
                    , data=SR20.22.9,REML=T)

SR.Sept.3reml<-lmer(log10(FLUX.mgCO2C.m2.h)~
                      RAINTRT * DISPERSION.FOURLEVELS * SPPNO
                    
                    #+ scale(log(MB.mg.gdw))
                    + scale(PROPTOTPLOT)  
                    + scale(lnRootBiomass)
                    + ScaledMoisture                              
                    + ScaledMoisture2
                    
                    + (1|PLOTNO) + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) 
                    , data=SR20.22.9,REML=T)

SR.Sept.4reml<-lmer(log10(FLUX.mgCO2C.m2.h)~
                      RAINTRT * DISPERSION.FOURLEVELS * SPPNO
                    
                    + scale(log(MB.mg.gdw))
                    + scale(PROPTOTPLOT)  
                    + scale(lnRootBiomass)
                    #+ ScaledMoisture                              
                    #+ ScaledMoisture2
                    
                    + (1|PLOTNO) + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) 
                    , data=SR20.22.9,REML=T)

SR.Sept.5reml<-lmer(log10(FLUX.mgCO2C.m2.h)~
                      RAINTRT * DISPERSION.FOURLEVELS * SPPNO
                    
                    + scale(log(MB.mg.gdw))
                    + scale(PROPTOTPLOT)  
                    #+ scale(lnRootBiomass)
                    + ScaledMoisture                              
                    + ScaledMoisture2
                    
                    + (1|PLOTNO) + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) 
                    , data=SR20.22.9,REML=T)


SR.Septreml<-lmer(log10(FLUX.mgCO2C.m2.h)~
                    RAINTRT * DISPERSION.FOURLEVELS * SPPNO
                  
                  + scale(log(MB.mg.gdw))
                  #+ scale(PROPTOTPLOT)    
                  #+ scale(log(DryWt.g.m2))
                  + ScaledMoisture                              
                  + ScaledMoisture2
                  
                  + (1|PLOTNO) + (1|SUBBLOCK) + (1|RAINTRT:SUBBLOCK) 
                  , data=SR20.22.9,REML=T)

SoilRespirationSept<-aictab(Cand.models <- list(
  "CVR + RB + MB + SM + SM2" = SR.Sept.1reml,
  "RB + MB + SM + SM2" = SR.Sept.2reml,
  "CVR + RB + SM + SM2" = SR.Sept.3reml,
  "CVR + RB + MB" = SR.Sept.4reml,
  "CVR + MB + SM + SM2" = SR.Sept.5reml,
  "MB + SM + SM2" = SR.Septreml))

SoilRespirationSept

r.squaredGLMM(SR.Sept.1reml)
r.squaredGLMM(SR.Sept.2reml)
r.squaredGLMM(SR.Sept.3reml)
r.squaredGLMM(SR.Sept.4reml)
r.squaredGLMM(SR.Sept.5reml)
r.squaredGLMM(SR.Septreml)

anova(SR.Sept<-lmer(log10(FLUX.mgCO2C.m2.h)~
                      RAINTRT                                
                    + DISPERSION.FOURLEVELS                   
                    + SPPNO
                    + RAINTRT:DISPERSION.FOURLEVELS 
                    + RAINTRT:SPPNO
                    + DISPERSION.FOURLEVELS:SPPNO
                    + RAINTRT:DISPERSION.FOURLEVELS:SPPNO
                    
                    + scale(log(MB.mg.gdw))
                    + scale(log(DryWt.g.m2))
                    + ScaledMoisture                              
                    + ScaledMoisture2
                    
                    + (1|PLOTNO) 
                    + (1|SUBBLOCK) 
                    + (1|RAINTRT:SUBBLOCK) 
                    , data=SR20.22.9),REML=T)

plot(SR.Sept)
anova(SR.Sept)
rand(SR.Sept)
r.squaredGLMM(SR.Sept)

