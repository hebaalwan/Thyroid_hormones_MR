#MR TSH UKBB phecodes

library(MendelianRandomization)
library(readxl)
library(dplyr)
#remotes::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)

###Loading data

merged_loci<- data.table::fread(".csv")

tsh.loci <- read_excel("")

tsh.loci<- select(tsh.loci, c("pos","SE"))
tsh.loci <- rename(tsh.loci, tsh.se=SE)

merged_loci <- merge(merged_loci, tsh.loci, by="pos")

detach()
attach(merged_loci)

#Thyroid cancer - date_193

res.presso.thycancer <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_193",
                        SdExposure = "tsh.se", SdOutcome = "SE.date_193", data=as.data.frame(merged_loci),
                        OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  

# Benign neoplasms of thyroid gland - date_226

res.presso.ben.neothy <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_226",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_226", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  


# Simple and unspecified goiter - date_240

res.presso.simplegoiter <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_240",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_240", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  

# Nontoxic uninodular goiter - date_241.1

res.presso.nontoxuni <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_241.1",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_241.1", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  

# Nontoxic multinodular goiter - date_241.2

res.presso.nontoxmultinod <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_241.2",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_241.2", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  

# Thyrotoxicosis with or without goiter - date_242

res.presso.thyrotox <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_242",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_242", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  

# Grave's disease - date_242.1

res.presso.graves <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_242.1",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_242.1", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  

# Toxic multinodular goiter - date_242.2

res.presso.multinodgoiter <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_242.2",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_242.2", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  


# Secondary hypothyroidism - date_244.1

res.presso.sechhypo <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_244.1",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_244.1", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  


# Acquired hypothyroidism - date_244.2

res.presso.acqhypo <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_244.2",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_244.2", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  

# Hypothyroidism NOS - date_244.4

res.presso.hypotnos <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_244.4",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_244.4", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  

# Other disorders of thyroid - date_246

res.presso.otherthyroid <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_246",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_246", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  

# Polydipsia - date_276.8

res.presso.polydipsia <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_276.8",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_276.8", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  


# Atrial fibrillation - date_427.21

res.presso.af <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_427.21",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_427.21", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  

# Hypercholesterolemia - date_272.11

res.presso.hyperchol <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_272.11",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_272.11", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  

# Pemphigus - date_695.22

res.presso.pemphigus <- mr_presso(BetaExposure = "beta.tsh", BetaOutcome = "beta.phecode.date_695.22",
                                  SdExposure = "tsh.se", SdOutcome = "SE.date_695.22", data=as.data.frame(merged_loci),
                                  OUTLIERtest = T, DISTORTIONtest = T, SignifThreshold = .05, NbDistribution = 5000, seed = 42)  





