#############################################################
#### PGS for thyroid function against UKBB phecodes       ####
####      ####################################################

rm(list = ls())
options(stringsAsFactors = F)
##load(".RData")

library(dplyr)
library(haven)

###############################################
####          load covariate data          ####
###############################################

## import basic covariate data created on linux
ukbb.cov       <- data.table::fread("/UKBB.covariates.csv", sep = ",", header=T, data.table=F)
## convert batch to binary variable
ukbb.cov$batch <- ifelse(ukbb.cov$batch < 0, 0, 1)

## import exlcusion list (white European unrelated samples)
ex.list        <- data.table::fread("/EXCLUDEFOR_White_Euro_Unrelateds_v1.samples", sep="\t", header=F)
## drop from the file
ukbb.cov       <- subset(ukbb.cov, !(n_eid_9905 %in% ex.list$V1))
## N = 352935


###############################################
####              import PGS              ####
###############################################


tsh.pgs <- data.table::fread("/tsh_loci_chr_ukbb_PGS.txt")

#-----------------------------------#
##--        add to the data      --##
#-----------------------------------#

## add to covariate data

ukbb.cov <- merge(ukbb.cov, tsh.pgs, by.x="n_eid_9905", by.y="IID")

###############################################
####            import phecodes            ####
###############################################

## load prepared phecode data
ukbb.phe <- data.table::fread("/UKBB.phecodes.collated.wide.covariates.20210922.txt", sep="\t", header=T, data.table=F)

## load label for phecodes, that is, mapping of column IDs to human readable names
lab.phe  <- read.delim("/Phecode.label.UKBB.20210922.txt",
                       sep="\t", header=T)

## keep only what is needed, drop some variables
ukbb.phe <- ukbb.phe[, c("f.eid", lab.phe$id)]

## import mapping table for 9905 IDs
map.tab  <- read_dta("/44448_9905_link_28Mar2019.dta")
## add to the phecode data
ukbb.phe <- merge(map.tab[, c("n_eid_9905", "n_eid_44448")], ukbb.phe, by.x="n_eid_44448", by.y="f.eid")

###Add scores
ukbb.phe <- merge(ukbb.cov, ukbb.phe, by.x="n_eid_9905", by.y="n_eid_9905")
rm(ukbb.cov); gc()


###############################################
####     run logistic regression models    ####
###############################################
## define adjustment set
adj    <- paste(c("age", "sex", "batch", "centre", paste0("pc", 1:10)), collapse = " + ")
pc <- paste0("pc", 1:10)


##Running logistic regression models for each phecode with PGS as the independent variable and showing NA for non-converged models
ukbb.phe$centre = factor(ukbb.phe$centre)

y=grep("^date", names(ukbb.phe), value=TRUE)   ## y is the list of variable names starting with ?date?
results = sapply(y, function(x){
  tryCatch({
    myfit=glm(formula= as.formula(paste0(x, "~ PGS + ", adj)), family = binomial, data= ukbb.phe ,na.action = na.omit)
    if(myfit$converged) {
      # estimates of pgs
      res1=summary(myfit)$coefficients["PGS",c(1,2,4)]
      # number of cases and controls used in the model
      ncase = nrow(myfit$model[myfit$model[,x]==1,])
      nctrl = nrow(myfit$model[myfit$model[,x]==0,])
      # Odds ratio and its calculate 95% CI
      OR=exp(res1[1])
      lowCI=exp(res1[1]-qnorm(0.975)*res1[2])
      upperCI=exp(res1[1]+qnorm(0.975)*res1[2])
      c(nctrl,ncase,res1,OR,lowCI,upperCI)
    }
    else {
      ncase = nrow(ukbb.phe[ukbb.phe[,x]==1,])
      nctrl = nrow(ukbb.phe[ukbb.phe[,x]==0,])
      c(nctrl,ncase,rep(NA,6))
    }
  }, error = function(e) {
    ncase = nrow(ukbb.phe[ukbb.phe[,x]==1,])
    nctrl = nrow(ukbb.phe[ukbb.phe[,x]==0,])
    c(nctrl,ncase,rep(NA,6))
  })
})
results=data.frame(t(results))
names(results)=c("n.controls", "n.cases", "beta", "SE","p", "OR", "lowCI", "upperCI")


#Removing phecodes with less than 200 cases

results       <- data.table::fread("/", sep = ",", header=T, data.table=F)

library(dplyr)

results <- subset(results, n.cases>200)

#Merge with labels file
results <- merge(lab.phe, results, by.x="id", by.y="V1")


#FDR
results <- results %>% mutate("adjusted.p.value" = p.adjust(p, method="fdr"))


