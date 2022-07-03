#Just loading some R packages.
library(foreign)
library(car)
library(plyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(stringr)
library(ggplot2)
library(jtools)
library(lmerTest)
library(ranger)
library(randomForest)
library(caret)
library(ROCR)
library(glmnet)
library(mice)
library(missRanger)
library(naniar)
library(grf)
library(data.table)
library(mltools)
library(bcf)
library(glmnet)
library(bartCause)
library(stats)

#Set random seed.
set.seed(123)

#Load the English subset of the data.
load("~/TALIS_DATA/TALIS_ENG.Rda")
df<-as.data.frame(df)

#Get final teacher weights.
mywgt<-as.numeric(as.character(df$TCHWGT))

#Remove some columns we don't need to use.
remove<-c(1:7, 352:460, 492:493)
df<-df[,-remove]


#Here we are converting categorical variables to continuous variables.
#Example: ("Not at all", "To some extent", "Quite a bit", "A lot")->(1, 2, 3, 4).
#This allows us to maintain the ordering of the categories and can be better than one-hot-encoding.
from=c("ENG1        ", "ENG2        ", 
       "Female", "Male", 
       "Not Reached", "Not administered", "Missing", "Omitted or invalid",
       "Logically not applicable",
       "Below <ISCED 2011 Level 3>", "<ISCED 2011 Level 3>", "<ISCED 2011 Level 4>",
       "<ISCED 2011 Level 5>", "<ISCED 2011 Level 6>", "<ISCED 2011 Level 7>",
       "<ISCED 2011 Level 8>", 
       "Yes", "No",
       "Not at all", "Somewhat", "Well", "Very well",
       "Not important at all", "Of low importance", "Of moderate importance",
       "Of high importance",
       "Full-time (more than 90% of full-time hours)", "Part-time (71-90% of full-time hours)",
       "Part-time (50-70% of full-time hours)", "Part-time (less than 50% of full-time hours)",
       "None", "Some", "Most", "All",
       "Checked", "Not checked",
       "No need at present", "Low level of need", "Moderate level of need", "High level of need",
       "Strongly disagree", "Disagree", "Agree", "Strongly agree",
       "Never", "Once a year or less", "2-4 times a year", "5-10 times a year", 
       "1-3 times a month", "Once a week or more",
       "Not at all", "To some extent", "Quite a bit", "A lot",
       "None", "1% to 10%", "11% to 30%", "31% to 60%", "More than 60%",
       "Never or almost never", "Occasionally", "Frequently", "Always")



to=c(1, 2,
     0, 1,
     NA, NA, NA, NA, NA,
     2, 3, 4, 5, 6, 7, 8,
     1, 0,
     1, 2, 3, 4,
     1, 2, 3, 4,
     90, 80, 60, 50,
     1, 2, 3, 4,
     1, 0,
     1, 2, 3, 4,
     1, 2, 3, 4,
     1, 2, 3, 4, 5, 6,
     1, 2, 3, 4,
     1, 5, 20, 50, 60,
     1, 2, 3, 4)


mymap<-function(x)
{
  x<-mapvalues(x,
               from=from,
               to=to)
  
  x<-as.numeric(as.character(x))
  
  return(x)
}


#Convert the categories to continuous variables.
df<-mutate_at(df, c(1:3, 6:37, 39:40, 45:71, 84:233, 239:291, 293:344), mymap)
df<-mutate_at(df, c(5, 41:44, 72:83, 235:238, 292, 345:375), as.character) 
df<-mutate_at(df, c(5, 41:44, 72:83, 235:238, 292, 345:375), as.numeric) 


#Use the missRanger package to impute missing data values.
df<-missRanger(df, pmm.k=3, num.trees=20)


#Select the final variables we want to keep and create a smaller subset of the data.
keep<-c(1:5, 30:45, 47, 234:238, 345:375)
dfsmall<-df[,keep]





#Create matrix of covariates.
X<-as.matrix(one_hot(as.data.table(dfsmall[,-c(6:12, 42, 43, 56)])))
cn<-colnames(X)
colnames(X)<-paste0("V", as.character(1:length(cn)))

#Create response variable: Desire to Move.
Y<-ifelse(df$TT3G53C>=3, 1, 0)

#Create treatment variable: 1 if teacher took part in Induction, 0 otherwise.
W<-ifelse(df$TT3G19A2==1|df$TT3G19B2==1, 1, 0)
W<-ifelse(is.na(W), 0, W)

#Run the bartCause model one initial time to get propensity score estimates.
tau.forest<-bartc(response=Y,
                  treatment=W,
                  confounders=X,
                  method.rsp="bart",
                  method.trt="bart",
                  estimand="ate",
                  commonSup.rule="none",
                  p.scoreAsCovariate=TRUE,
                  weights=mywgt,
                  ntree=350,
                  nskip=4000,
                  ndpost=4000,
                  n.chains=1,
                  keepTrees=TRUE,
                  args.trt=c(k=300))


#The propensity score estimates.
W.hat<-fitted(tau.forest, "p.score")



#Load the TALIS data again because we removed the weights earlier but now we need them for the BRR procedure.
load("~/TALIS_DATA/TALIS_ENG.Rda")
df<-as.data.frame(df)



#Create a list for storing the results obtained from using each different weight.
av_nums<-c()

#Loop through each of the 100 weights.
for(i in 1:100)
{
  
  tau.forest<-bartc(response=Y,
                  treatment=W,
                  confounders=X,
                  method.rsp="bart",
                  #Here we use the propensity scores estimated earlier.
                  method.trt=W.hat,
                  estimand="ate",
                  commonSup.rule="none",
                  p.scoreAsCovariate=TRUE,
                  #Select the correct index for the weights below.
                  weights=as.numeric(as.character(df[,359+i])),
                  ntree=350,
                  nskip=5000,
                  ndpost=5000,
                  n.threads=6)

  #Get individual treatment effect estimates.
  tau.hat<-fitted(tau.forest, type="icate")
  
  #Append the resulting ATE estimate to the list.
  av_nums<-c(av_nums, weighted.mean(tau.hat, as.numeric(as.character(df[,359+i]))))
}

#Calculate the standard error.
av<-mean(av_nums)
vfay1<-sum((av_nums-av)^2)
vfay2<-vfay1/(100*(1-0.5)^2)
se<-vfay2^0.5

#Print the ATE and standard error.
paste0("ATE: ", av, " ", "SE: ", se)



