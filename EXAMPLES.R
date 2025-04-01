

source("./R/FUNCTIONS_EXAMPLE.R") #load functions written for this workshop

library(ShinyItemAnalysis) # Grant example data
library(reshape2)          # to transform data from long to wide etc.
library(rjags)             # Fit models in JAGS
library(irr)               # Compute ICC
library(lme4)              #Frequentist multilevel probit model


#######################################
#
#GRANT REVIEW DATA ANALYSIS
#
#######################################



##############################
# Creating the dataset
##############################

data(AIBS)
data.long<-AIBS[AIBS[,2]==2015,c(1,2,13,25)]
data.long[,4]<-factor(data.long[,4])

# Wide dataset
###############
data.wide<-data.frame(dcast(data.long,ID~RevCode,value.var=c("Score")))

#Wide dataset - binary
#######################
databin.wide<-matrix(NA,ncol=4,nrow=nrow(data.wide))
for (i in 1:3){databin.wide[,i+1]<-ifelse(data.wide[,i+1]<=2,1,0)}
databin.wide[,1]<-seq(1,29)
N<-nrow(databin.wide) #number of subjects
table(databin.wide[,c(2,3,4)])/(29*3)

################################################
#kappa coefficient - chance 3 because ANOVA1
################################################

kappaw(2-databin.wide[,c(2,3,4)],weight="binary",chance=3)

kappa.bayes.bin(2-databin.wide[,c(2,3,4)])

##############################
# Multilevel probit model
##############################

#Transform in long dataset
############################
n<-nrow(databin.wide[,c(2,3,4)])
nrat<-ncol(databin.wide[,c(2,3,4)])

y<-c(databin.wide[,c(2,3,4)]) 
sub<-c(rep(seq(1,n),nrat))
dat<-data.frame(cbind(y,sub))

#Bayesian
############
  data.jags<-list(n=nrat*n,ns=n,y=y,sub=sub)
  parameters<-c("ICC","phi","sigma2_S")
  jags.inits<-list(beta0=0.0,sigma_S=4)
  obj <- jags.model("probit5.txt",data=data.jags,n.chains=3,inits=jags.inits,quiet=TRUE)
  jags.samples(obj, parameters, n.iter=2000, progress.bar="none");
  model.s<-coda.samples(obj, parameters, n.iter=5000, progress.bar="none")
  summary(model.s)
 
# Frequentist
#############
l_res<-glmer(y~1+(1|sub),data=dat,family='binomial'(link = "probit"),control = glmerControl(optimizer = "bobyqa"),nAGQ = 15)

probICC<-prob_icc1(l_res,'sub')
    
b.ICC<-bootp_icc1(l_res, subject, 1000, seed=2025)
    
list("ICC.LAT"=probICC[1],"LOW.LAT"=b.ICC[1,1],"HIGH.LAT"=b.ICC[3,1], "ICC"=probICC[2],"LOW"=b.ICC[1,2],"HIGH"=b.ICC[3,2])
 


#######################################
#
#ENGLISH CLASS EXAMPLE
#
#######################################   

data.t<-read.csv("./LANG2.csv",sep=";",header=TRUE)

#Randomly select one teacher
set.seed(2025)
sample(1:3, 1)

#Data for teacher 1
teacher1<-data.t[data.t[,1]==1,][,-c(1,2)]
teacher2<-data.t[data.t[,1]==2,][,-c(1,2)]
teacher3<-data.t[data.t[,1]==3,][,-c(1,2)]

#Transform to binary data
teacher1<-ifelse(teacher1==3|teacher1==4,1,0)
teacher2<-ifelse(teacher2==3|teacher2==4,1,0)
teacher3<-ifelse(teacher3==3|teacher3==4,1,0)

#Proportion of cognitive tasks 
colMeans(teacher1)
  
#kappa coefficient - chance 2 because reproducibility study
###############################################################

kappaw(teacher1+1,"binary",2)
kappa.bayes.bin(data.frame(teacher1+1))


#Multilevel probit model - reproducibility study
###############################################################

#Transform the dataset into long dataset
########################################
n<-nrow(teacher1)
nrat<-ncol(teacher1)

y<-c(teacher1) 
sub<-c(rep(seq(1,n),nrat))
met<-sort(rep(1:nrat,n))
dat<-data.frame(cbind(y,sub,met))

# Bayesian analysis - Cauchy prior - reproducibility
#####################################################

data.jags<-list(n=nrat*n,ns=n,nr=nrat,y=y,sub=sub,met=met)
parameters<-c("ICC","beta0","sigma2_S","sigma2_R","phi")
jags.inits<-list(beta0=0.5,sigma_S=2,sigma_R=0.1)
obj <- jags.model("probit6.txt",data=data.jags,n.chains=3,inits=jags.inits,quiet=TRUE)
jags.samples(obj, parameters, n.iter=2000, progress.bar="none");
model.s<-coda.samples(obj, parameters, n.iter=5000, progress.bar="none")
summary(model.s)

#Multilevel probit model - reproducibility study -Frequentist
###############################################################

l_res<-glmer(y~1+(1|sub)+(1|met),data=dat,family='binomial'(link = "probit"))

probICC<-prob_icc2(l_res,'sub','met')
  
b.ICC<-bootp_icc2(l_res, subject, B=1000, seed=2025)

#list("ICC.LAT"=b.ICC[2,1],"LOW.LAT"=b.ICC[1,1],"HIGH.LAT"=b.ICC[3,1],"ICC"=b.ICC[2,2],"LOW"=b.ICC[1,2],"HIGH"=b.ICC[3,2])



#######################################
#
#FETAL HEART RATE DATA DATA ANALYSIS
#
#######################################

#Read the data
##################################
HR<-read.csv("./HR.csv")

#Make the binary variable
##################################
HR_01<-ifelse(HR_NA>110 &HR_NA<150,0,1)

#Proportion of normal/abnormal
##################################
colMeans(HR_01)

#Absolute ICC two-way ANOVA (Frequentist)
##################################
icc(HR,model="twoway",type="agreement",unit="single",conf.level=0.95)


# Probit - Bayesian (Cauchy prior)
##################################

#Transform the data in long dataset
  n<-nrow(HR_01)
  nrat<-ncol(HR_01)
  
  y<-c(HR_01) 
  sub<-c(rep(seq(1,n),nrat))
  met<-sort(rep(1:nrat,n))

#Fit the multilevel probit model
  data.jags<-list(n=nrat*n,ns=n,nr=nrat,y=y,sub=sub,met=met)
  parameters<-c("ICC","beta0","sigma2_S","sigma2_R","phi")
  jags.inits<-list(beta0=0.5,sigma_S=2,sigma_R=0.1)
  obj <- jags.model("probit6.txt",data=data.jags,n.chains=3,inits=jags.inits,quiet=TRUE)
  jags.samples(obj, parameters, n.iter=2000, progress.bar="none");
  model.s<-coda.samples(obj, parameters, n.iter=5000, progress.bar="none")
  summary(model.s)
  






