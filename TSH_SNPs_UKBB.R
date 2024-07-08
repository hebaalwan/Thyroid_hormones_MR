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


tsh.snps <- read.table(gzfile("/tsh_loci_chr_ukbb_bysnp.dosage.gz"), header=T)
snp.info <- read.table("/tsh_loci_chr_ukbb_bysnp.stats", header=T)

#tsh.snps = subset(tsh.snps, select = -c(rs12390237))

#-----------------------------------#
##--        add to the data      --##
#-----------------------------------#

## add to covariate data

ukbb.cov <- merge(ukbb.cov, tsh.snps, by.x="n_eid_9905", by.y="IID")

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

ukbb.phe <- ukbb.phe[, c("f.eid", "date_193", "date_226", "date_240", "date_241.1", "date_241.2", "date_242", "date_242.1", "date_242.2", "date_244.1","date_244.2", "date_244.4", "date_246", "date_272.11", "date_276.8", "date_427.21", "date_695.22"
)]

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


adj    <- paste(c("age", "sex", "batch", "centre", paste0("pc", 1:10)), collapse = " + ")


##Running logistic regression models for each phecode with PGS as the independent variable and showing NA for non-converged models
ukbb.phe$centre = factor(ukbb.phe$centre)

y=grep("^date", names(ukbb.phe), value=TRUE)   ## y is the list of variable names starting with ?date?

## make it faster by doing it in parallel
require(doMC)
registerDoMC(12)

marker.list=names(tsh.snps)[2:dim(tsh.snps)[2]]  # list of SNPs
myres = lapply(marker.list, function(m) {
  results = sapply(y, function(x){
    tryCatch({
      myfit=glm(formula= as.formula(paste0(x, "~ ", m, " + ", adj)), family = binomial, data= ukbb.phe ,na.action = na.omit)
      if(myfit$converged) {
        # estimates of pgs
        res1=summary(myfit)$coefficients[2,c(1,2,4)]
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
  results %>% mutate("adjusted p-value" = p.adjust(p, method="fdr"))
})


myres         <- `names<-`(myres, marker.list)  # Rename the component of the list
myres         <- do.call(rbind, myres)
## add SNP exposure
myres$snp     <- sapply(row.names(myres), function(x) strsplit(x, "\\.")[[1]][1])
myres$phecode <- sapply(row.names(myres), function(x) paste0("date",strsplit(x, "\\.date")[[1]][2]))

## add allele information
myres         <- merge(myres, snp.info, by.x="snp", by.y="rsid")

## add TSH loci
tsh.loci      <- read_dta("/tsh_loci.dta")

## convert chromosome to numeric
myres$chr = ifelse(myres$chr=="X", 23, myres$chr)
myres$chr     <- as.numeric(myres$chr)

## add to wide
myres.wide    <- merge(tsh.loci, myres, by.x=c("rsid", "chr", "pos"), by.y = c("snp", "chr", "pos"), suffixes = c(".tsh", ".phecode"))

## align alleles
myres.wide$beta.phecode <- ifelse(myres.wide$effect_allele == myres$alleleB, myres.wide$beta.phecode, -myres.wide$beta.phecode)
## drop alleles no longer needed
myres.wide$alleleA <- myres.wide$alleleB <- NULL

## make all effects coded in positive direction for exposure (TSH)
myres.wide$beta.phecode      <- sign(myres.wide$beta.tsh)*myres.wide$beta.phecode
## generate aligning allele definition
myres.wide$effect_allele_new <- ifelse(sign(myres.wide$beta.tsh) == 1, myres.wide$effect_allele, myres.wide$other_allele)
myres.wide$other_allele_new  <- ifelse(sign(myres.wide$beta.tsh) != 1, myres.wide$effect_allele, myres.wide$other_allele)
## now convert all TSH association to be positive
myres.wide$beta.tsh          <- abs(myres.wide$beta.tsh)
## N.B.: alleles information might no longer be valid


## change to wide format
myres.wide    <- reshape(myres.wide, idvar=c("rsid", "chr", "pos", "HWEpval", "FreqB", "info", "other_allele", "effect_allele", "beta.tsh"),
                         timevar = "phecode", direction = "wide")

