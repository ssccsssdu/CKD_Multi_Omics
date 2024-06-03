
library(TwoSampleMR)
library(data.table)
library(openxlsx)
library(dplyr)
library(fdrtool)
library(meta)
library(ggforce)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(Cairo)
library(ggpubr)
library(ggvenn)
library(Cairo)
library(corrplot)
library(coloc)

###############################################PWAS analysis####################################################
################################################################################################################

#####################################read and format CKD GWAS data##############################################

ckd <- fread("CKD_overall_EA_JW_20180223_nstud23.dbgap.txt")

ckd$p <- ckd$`P-value`

ckd <- format_data(
  ckd,
  type='outcome',
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col ="Allele1",
  other_allele_col = "Allele2",
  pval_col = "p",
  eaf_col = "Freq1",
  samplesize_col = "n_total_sum",
  units_col='cis', 
  min_pval=0
)

ckd$outcome <- 'CKD'
ckd$id.outcome <- 'CKD'
ckd$units.outcome <- 'log odds'
ckd$prevalence.outcome <- 41395/439303
ckd$ncase.outcome <- 41395
ckd$ncontrol.outcome <- 439303-41395
ckd$samplesize.outcome <- 439303

##########################################pqtl 1: read and format deCODE pqtl#################################################

pqtl1 <- read.xlsx('deCODE_pqtl.xlsx')

pqtl1$se <- abs(pqtl1$b/qnorm((10^-(pqtl1$log10p))/2))
pqtl1$p <- 10^-(pqtl1$log10p)
pqtl1$exposure <- pqtl1$target
pqtl1$maf <- pqtl1$maf/100
pqtl1 <- format_data(
  pqtl1,
  type='exposure',
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col ="A1",
  other_allele_col = "A0",
  pval_col = "p",
  eaf_col = "eaf",
  phenotype_col='target',
  min_pval=0,
  units_col='cis',
  chr_col = "chr",
  pos_col = "pos"
)

pqtl1$type <- pqtl1$units.exposure
pqtl1$units.exposure <- 'SD'
pqtl1$samplesize.exposure <- 35597
pqtl1 <- pqtl1[-grep(',',pqtl1$SNP),]  #exclude snp without snpid
pqtl1 <- pqtl1[-grep('chr',pqtl1$SNP),] #exclude snp without snpid
pqtl1 <- pqtl1[!is.na(pqtl1$effect_allele.exposure)&!is.na(pqtl1$other_allele.exposure),]

#######################################pqtl 2: read and format UKB_PPP pqtl######################################################

pqtl2 <- read.xlsx('UKB_PPP_pqtl.xlsx')
pqtl2$p <- 10^-(pqtl2$log10p)
names(pqtl2)
pqtl2$eaf <- ifelse(pqtl2$eaf>0.5,1-pqtl2$eaf,pqtl2$eaf)

pqtl2 <- format_data(
  pqtl2,
  type='exposure',
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col ="A1",
  other_allele_col = "A0",
  pval_col = "p",
  eaf_col = "eaf",
  phenotype_col='target',
  units_col='cis',
  min_pval=0,
  chr_col = "chr",
  pos_col = "pos"
)

pqtl2$type <- pqtl2$units.exposure
pqtl2$units.exposure <- 'SD'
pqtl2$samplesize.exposure <- 54219

#######################################pqtl 3: read and format Fenland pqtl ######################################################

pqtl3 <- read.xlsx('Fenland_pqtl.xlsx')
names(pqtl3)

pqtl3 <- format_data(
  pqtl3,
  type='exposure',
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col ="A1",
  other_allele_col = "A0",
  pval_col = "p",
  eaf_col = "eaf",
  phenotype_col='target',
  units_col='cis',
  min_pval=0,
  chr_col = "chr",
  pos_col = "pos"
)

pqtl3$type <- pqtl3$units.exposure
pqtl3$units.exposure <- 'SD'
pqtl3$samplesize.exposure <- 10708

#####################################harmonise pqtl data with CKD##############################################

dat1 <- harmonise_data(pqtl1,ckd)
dat1 <- dat1[dat1$mr_keep,]
dat1 <- dat1[!is.na(dat1$eaf.outcome),]
dat1 <- steiger_filtering(dat1)
dat1 <- dat1[dat1$steiger_dir,]
dat1_cis <- dat1[dat1$type=='cis',]

#######################

dat2 <- harmonise_data(pqtl2,ckd)
dat2 <- dat2[dat2$mr_keep,]
dat2 <- dat2[!is.na(dat2$eaf.outcome),]
dat2 <- steiger_filtering(dat2)
dat2 <- dat2[dat2$steiger_dir,]
dat2_cis <- dat2[dat2$type=='cis',]

#######################

dat3 <- harmonise_data(pqtl3,ckd)
dat3 <- dat3[dat3$mr_keep,]
dat2 <- dat2[!is.na(dat2$eaf.outcome),]
dat3 <- steiger_filtering(dat3)
dat3 <- dat3[dat3$steiger_dir,]
dat3_cis <- dat3[dat3$type=='cis',]

######################################principal PWAS analyis##############################################################

res1 <- mr(dat1_cis,method_list = c('mr_wald_ratio','mr_ivw'))
res1$exposure <- gsub(' (cis)','',res1$exposure,fixed = T)

res2 <- mr(dat2_cis,method_list = c('mr_wald_ratio','mr_ivw'))
res2$exposure <- gsub(' (cis)','',res2$exposure,fixed = T)

res3 <- mr(dat3_cis,method_list = c('mr_wald_ratio','mr_ivw'))
res3$exposure <- gsub(' (cis)','',res3$exposure,fixed = T)

###########select proteins passed FDR correction

res1$dat <- 'Iceland'
res2$dat <- 'UK Biobank'
res3$dat <- 'Fenland'

res <- rbind(res1,res2,res3)
fdr <- fdrtool(res$pval,statistic="pvalue",plot=F)
res$qval <- fdr$qval 
res.fdr <- res[res$qval<0.05,]

######output cis data

dat1_cis$exposure <- gsub(' (cis)','',dat1_cis$exposure,fixed = T)
dat1_cis2 <- dat1_cis[dat1_cis$exposure%in%res.fdr$exposure,]
dat1_cis2$data <- 'Iceland'
dat2_cis$exposure <- gsub(' (cis)','',dat2_cis$exposure,fixed = T)
dat2_cis2 <- dat1_cis[dat2_cis$exposure%in%res.fdr$exposure,]
dat2_cis2$data <- 'UK Biobank'
dat3_cis$exposure <- gsub(' (cis)','',dat3_cis$exposure,fixed = T)
dat3_cis2 <- dat1_cis[dat3_cis$exposure%in%res.fdr$exposure,]
dat3_cis2$data <- 'Fenland'

dat_cis <- bind_rows(dat1_cis2,dat2_cis2,dat3_cis2)
#write.xlsx(dat_cis,'dat_cis.xlsx')

##############combined effects by Meta analysis

res1.sig <- res1[res1$exposure%in%res.fdr$exposure,]
res1.sig$group <- 'cis'
res1.sig$data <- 'Iceland'

res2.sig <- res2[res2$exposure%in%res.fdr$exposure,]
res2.sig$group <- 'cis'
res2.sig$data <- 'UK Biobank'

res3.sig <- res3[res3$exposure%in%res.fdr$exposure,]
res3.sig$group <- 'cis'
res3.sig$data <- 'Fenland'

res.sig <- rbind(res1.sig,res2.sig,res3.sig)
res.sig <- res.sig%>%group_by(exposure)%>%mutate(num=length(exposure))%>%ungroup()
res.sig2 <- res.sig[res.sig$num>1,] #protein in more than 1 dataset

res.combine <- NULL

for(i in unique(res.sig2$exposure)){
  dt <- res.sig2[res.sig2$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt$source <- paste(dt$data,collapse = ',')
  dt <- dt[1,]
  dt$data <- 'Combined'
  res.combine <- rbind(res.combine,dt)
}

res.sig <- bind_rows(res.sig,res.combine)
res.sig <- generate_odds_ratios(res.sig)
res.sig$CI <- paste0(format(round(res.sig$or,2),nsmall = 2),'(',format(round(res.sig$or_lci95,2),nsmall = 2),'-',format(round(res.sig$or_uci95,2),nsmall = 2),')')

##################results with cis+trans pqtl

dat1$exposure <- gsub(' (cis)','',dat1$exposure,fixed = T)
dat1$exposure <- gsub(' (trans)','',dat1$exposure,fixed = T)
dat1a <- dat1[dat1$exposure%in%res1.sig$exposure,]
dat1a$id.exposure <- dat1a$exposure
res1a <- mr(dat1a,method_list = c('mr_wald_ratio','mr_ivw'))
res1a$group <- 'cis_trans'
res1a$data <- 'Iceland'

dat2$exposure <- gsub(' (cis)','',dat2$exposure,fixed = T)
dat2$exposure <- gsub(' (trans)','',dat2$exposure,fixed = T)
dat2a <- dat2[dat2$exposure%in%res2.sig$exposure,]
dat2a$id.exposure <- dat2a$exposure
res2a <- mr(dat2a,method_list = c('mr_wald_ratio','mr_ivw'))
res2a$group <- 'cis_trans'
res2a$data <- 'UK Biobank'

dat3$exposure <- gsub(' (cis)','',dat3$exposure,fixed = T)
dat3$exposure <- gsub(' (trans)','',dat3$exposure,fixed = T)
dat3a <- dat3[dat3$exposure%in%res3.sig$exposure,]
dat3a$id.exposure <- dat3a$exposure
res3a <- mr(dat3a,method_list = c('mr_wald_ratio','mr_ivw'))
res3a$group <- 'cis_trans'
res3a$data <- 'Fenland'

res.a <- rbind(res1a,res2a,res3a)
res.a <- res.a%>%group_by(exposure)%>%mutate(num=length(exposure))%>%ungroup()
res.a.2 <- res.a[res.a$num>1,]

res.a.combine <- NULL

for(i in unique(res.a.2$exposure)){
  dt <- res.a.2[res.a.2$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt$source <- paste(dt$data,collapse = ',')
  dt <- dt[1,]
  dt$data <- 'Combined'
  res.a.combine <- rbind(res.a.combine,dt)
}

res.a <- bind_rows(res.a,res.a.combine)
res.a <- generate_odds_ratios(res.a)
res.a$CI <- paste0(format(round(res.a$or,2),nsmall = 2),'(',format(round(res.a$or_lci95,2),nsmall = 2),'-',format(round(res.a$or_uci95,2),nsmall = 2),')')

########################MR Egger and heterogeneity test 

dat1a_cis <- dat1a[dat1a$type=='cis',]
dat2a_cis <- dat2a[dat2a$type=='cis',]
dat3a_cis <- dat3a[dat3a$type=='cis',]

egger1 <- mr_pleiotropy_test(dat1a_cis)
egger1$dat <- 'Iceland'
egger2 <- mr_pleiotropy_test(dat2a_cis)
egger2$dat <- 'UK Biobank'
egger3 <- mr_pleiotropy_test(dat3a_cis)
egger3$dat <- 'Fenland'
egger <- rbind(egger1,egger2,egger3)
write.xlsx(egger,'egger.xlsx')

hetero1 <- mr_heterogeneity(dat1a_cis)
hetero1$dat <- 'Iceland'
hetero2 <- mr_heterogeneity(dat2a_cis)
hetero2$dat <- 'UK Biobank'
hetero3 <- mr_heterogeneity(dat3a_cis)
hetero3$dat <- 'Fenland'
hetero <- rbind(hetero1,hetero2,hetero3)

write.xlsx(hetero,'hetero.xlsx')

###############################################TWAS results####################################################
################################################################################################################

#mapping all coding genes
pqtl1_gene <- read.xlsx('deCODE_pqtl.xlsx')
pqtl1_gene <- unique(pqtl1_gene[pqtl1_gene$target%in%res.fdr$exposure,'gene'])

pqtl2_gene <- read.xlsx('UKB_PPP_pqtl.xlsx')
pqtl2_gene <- unique(pqtl2_gene[pqtl2_gene$target%in%res.fdr$exposure,'gene'])

pqtl3_gene <- read.xlsx('Fenland_pqtl.xlsx')
pqtl3_gene <- unique(pqtl3_gene[pqtl3_gene$Target%in%res.fdr$exposure,'gene'])

all_gene <- c(pqtl1_gene,pqtl2_gene,pqtl3_gene)
all_gene <- unique(unlist(strsplit(all_gene,'|',fixed = T)))

gene_list <- read.table("2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",header = T)
gene_list <- unique(gene_list[,c('Gene','GeneSymbol')]) # acquire the GeneSymbol

###################eQTLGen results

gene1 <- read.table("ckd.smr",header=T) #input results of SMR analysis by smr software
gene1 <- merge(gene1,gene_list,by='Gene')
gene1 <- gene1[gene1$GeneSymbol%in%all_gene,]
gene1$Gene <- gene1$GeneSymbol

#######################brain results

gene2 <- read.table("ckd_brain.smr",header=T)
gene2 <- gene2[gene2$Gene%in%all_gene,]

#######################GTEx resullts

file <- list.files("GTEx",full.names = T)
file2 <- list.files("GTEx")

file <- file[grep('.smr',file,fixed=T)]
file2 <- file2[grep('.smr',file2,fixed=T)]
file2 <- gsub('.smr','',file2,fixed=T)

gene3 <- NULL
for(i in 1:49){
  dt <- fread(file[i])
  dt$tissue <- file2[i]
  gene3 <- rbind(gene3,dt)
}

gene3 <- gene3[gene3$Gene%in%all_gene,]

######################CAGE and westra results

gene4 <- read.table("ckd_CAGE.smr",header=T)
gene4 <- gene4[gene4$Gene%in%all_gene,]
gene4 <- gene4[gene4$probeID!='ILMN_1676528',] #only show one probe

gene5 <- read.table("ckd_westra.smr",header=T)
gene5 <- gene5[gene5$Gene%in%all_gene,]

#######################combined all results

gene1$tissue <- 'eQTLGen(Blood)'
gene2$tissue <- 'PsychENCODE(Brain)'
gene3$tissue <- paste0('GTEx(',gene3$tissue,")")
gene4$tissue <- 'CAGE(Blood)'
gene5$tissue <- 'Westra(Blood)'

gene12345 <- dplyr::bind_rows(gene1,gene4,gene5,gene2,gene3,)

#########################################examples of sensitivity and replication analysis#######################
####################################################################################################
 
##########CKD2  CKDGen(trans-ancestry,2019)

ckd2 <- fread("ckd_overall_ALL_JW_20180223_nstud30.dbgap.txt")

ckd2$p <- ckd2$`P-value`
ckd2$A1 <- toupper(ckd2$Allele1)
ckd2$A0 <- toupper(ckd2$Allele2)
head(ckd2)

ckd2 <- format_data(
  ckd2,
  type='outcome',
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col ="Allele1",
  other_allele_col = "Allele2",
  pval_col = "p",
  eaf_col = "Freq1",
  samplesize_col = "n_total_sum",
  units_col='cis',
  min_pval=0
)

ckd2$outcome <- 'ckd2'
ckd2$id.outcome <- 'ckd2'
ckd2$units.outcome <- 'log odds'
ckd2$prevalence.outcome <- 64164/625219
ckd2$ncase.outcome <- 64164
ckd2$ncontrol.outcome <- 625219-64164
ckd2$samplesize.outcome <- 625219

dat1.2 <- harmonise_data(pqtl1,ckd2)
dat1.2 <- dat1.2[dat1.2$type=='cis',]
dat1.2$exposure <- gsub(' (cis)','',dat1.2$exposure,fixed = T)
dat1.2 <- dat1.2[dat1.2$exposure%in%res.fdr$exposure,]
dat1.2 <- dat1.2[dat1.2$mr_keep,]
dat1.2 <- dat1.2[!is.na(dat1.2$eaf.outcome),]
dat1.2 <- steiger_filtering(dat1.2)
dat1.2 <- dat1.2[dat1.2$steiger_dir,]
res1.2 <- mr(dat1.2,method_list = c('mr_wald_ratio','mr_ivw'))
res1.2$exp <- 'Iceland'
res1.2$out <- 'CKDGen(trans-ancestry,2019)'

dat2.2 <- harmonise_data(pqtl2,ckd2)
dat2.2 <- dat2.2[dat2.2$type=='cis',]
dat2.2$exposure <- gsub(' (cis)','',dat2.2$exposure,fixed = T)
dat2.2 <- dat2.2[dat2.2$exposure%in%res.fdr$exposure,]
dat2.2 <- dat2.2[dat2.2$mr_keep,]
dat2.2 <- dat2.2[!is.na(dat2.2$eaf.outcome),]
dat2.2 <- steiger_filtering(dat2.2)
dat2.2 <- dat2.2[dat2.2$steiger_dir,]
res2.2 <- mr(dat2.2,method_list = c('mr_wald_ratio','mr_ivw'))
res2.2$exp <- 'UK Biobank'
res2.2$out <- 'CKDGen(trans-ancestry,2019)'

dat3.2 <- harmonise_data(pqtl3,ckd2)
dat3.2 <- dat3.2[dat3.2$type=='cis',]
dat3.2$exposure <- gsub(' (cis)','',dat3.2$exposure,fixed = T)
dat3.2 <- dat3.2[dat3.2$exposure%in%res.fdr$exposure,]
dat3.2 <- dat3.2[dat3.2$mr_keep,]
dat3.2 <- dat3.2[!is.na(dat3.2$eaf.outcome),]
dat3.2 <- steiger_filtering(dat3.2)
dat3.2 <- dat3.2[dat3.2$steiger_dir,]
res3.2 <- mr(dat3.2,method_list = c('mr_wald_ratio','mr_ivw'))
res3.2$exp <- 'Fenland'
res3.2$out <- 'CKDGen(trans-ancestry,2019)'

res_ckd1 <- rbind(res1.2,res2.2,res3.2)

#########################CKD3-4 

pqtl1a <- pqtl1
pqtl1a <- pqtl1a[pqtl1a$type=='cis',]
pqtl1a$exposure <- gsub(' (cis)','',pqtl1a$exposure,fixed = T)
pqtl1a <- pqtl1a[pqtl1a$exposure%in%res.fdr$exposure,]
ckd3 <- extract_outcome_data(snps=pqtl1a$SNP,outcomes=c('ebi-a-GCST003374','ebi-a-GCST90018822'))
dat1.3 <- harmonise_data(pqtl1a,ckd3)
dat1.3 <- dat1.3[dat1.3$mr_keep,]
dat1.3 <- dat1.3[!is.na(dat1.3$eaf.outcome),]
dat1.3$units.outcome <- 'log odds'
dat1.3$ncase.outcome <- NA
dat1.3$ncontrol.outcome <- NA
dat1.3$ncase.outcome <- ifelse(dat1.3$id.outcome=='ebi-a-GCST003374',12385,dat1.3$ncase.outcome)
dat1.3$ncontrol.outcome <- ifelse(dat1.3$id.outcome=='ebi-a-GCST003374',104780,dat1.3$ncontrol.outcome)
dat1.3$ncase.outcome <- ifelse(dat1.3$id.outcome=='ebi-a-GCST90018822',8287,dat1.3$ncase.outcome)
dat1.3$ncontrol.outcome <- ifelse(dat1.3$id.outcome=='ebi-a-GCST90018822',474571,dat1.3$ncontrol.outcome)
dat1.3$samplesize.outcome <- dat1.3$ncase.outcome+dat1.3$ncontrol.outcome
dat1.3$prevalence.outcome <- dat1.3$ncase.outcome/dat1.3$samplesize.outcome
dat1.3 <- steiger_filtering(dat1.3)
dat1.3 <- dat1.3[dat1.3$steiger_dir,]
res1.3 <- mr(dat1.3,method_list = c('mr_wald_ratio','mr_ivw'))
res1.3$exp <- 'Iceland'

##
pqtl2a <- pqtl2
pqtl2a <- pqtl2a[pqtl2a$type=='cis',]
pqtl2a$exposure <- gsub(' (cis)','',pqtl2a$exposure,fixed = T)
pqtl2a <- pqtl2a[pqtl2a$exposure%in%res.fdr$exposure,]
ckd3 <- extract_outcome_data(snps=pqtl2a$SNP,outcomes=c('ebi-a-GCST003374','ebi-a-GCST90018822'))
dat2.3 <- harmonise_data(pqtl2a,ckd3)
dat2.3 <- dat2.3[dat2.3$mr_keep,]
dat2.3 <- dat2.3[!is.na(dat2.3$eaf.outcome),]
dat2.3$units.outcome <- 'log odds'
dat2.3$ncase.outcome <- NA
dat2.3$ncontrol.outcome <- NA
dat2.3$ncase.outcome <- ifelse(dat2.3$id.outcome=='ebi-a-GCST003374',12385,dat2.3$ncase.outcome)
dat2.3$ncontrol.outcome <- ifelse(dat2.3$id.outcome=='ebi-a-GCST003374',104780,dat2.3$ncontrol.outcome)
dat2.3$ncase.outcome <- ifelse(dat2.3$id.outcome=='ebi-a-GCST90018822',8287,dat2.3$ncase.outcome)
dat2.3$ncontrol.outcome <- ifelse(dat2.3$id.outcome=='ebi-a-GCST90018822',474571,dat2.3$ncontrol.outcome)
dat2.3$samplesize.outcome <- dat2.3$ncase.outcome+dat2.3$ncontrol.outcome
dat2.3$prevalence.outcome <- dat2.3$ncase.outcome/dat2.3$samplesize.outcome
dat2.3 <- steiger_filtering(dat2.3)
dat2.3 <- dat2.3[dat2.3$steiger_dir,]
res2.3 <- mr(dat2.3,method_list = c('mr_wald_ratio','mr_ivw'))
res2.3$exp <- 'UK Biobank'

##
pqtl3a <- pqtl3
pqtl3a <- pqtl3a[pqtl3a$type=='cis',]
pqtl3a$exposure <- gsub(' (cis)','',pqtl3a$exposure,fixed = T)
pqtl3a <- pqtl3a[pqtl3a$exposure%in%res.fdr$exposure,]
ckd3 <- extract_outcome_data(snps=pqtl3a$SNP,outcomes=c('ebi-a-GCST003374','ebi-a-GCST90018822'))

dat3.3 <- harmonise_data(pqtl3a,ckd3)
dat3.3 <- dat3.3[dat3.3$mr_keep,]
dat3.3 <- dat3.3[!is.na(dat3.3$eaf.outcome),]
dat3.3$units.outcome <- 'log odds'
dat3.3$ncase.outcome <- NA
dat3.3$ncontrol.outcome <- NA
dat3.3$ncase.outcome <- ifelse(dat3.3$id.outcome=='ebi-a-GCST003374',12385,dat3.3$ncase.outcome)
dat3.3$ncontrol.outcome <- ifelse(dat3.3$id.outcome=='ebi-a-GCST003374',104780,dat3.3$ncontrol.outcome)
dat3.3$ncase.outcome <- ifelse(dat3.3$id.outcome=='ebi-a-GCST90018822',8287,dat3.3$ncase.outcome)
dat3.3$ncontrol.outcome <- ifelse(dat3.3$id.outcome=='ebi-a-GCST90018822',474571,dat3.3$ncontrol.outcome)
dat3.3$samplesize.outcome <- dat3.3$ncase.outcome+dat3.3$ncontrol.outcome
dat3.3$prevalence.outcome <- dat3.3$ncase.outcome/dat3.3$samplesize.outcome
dat3.3 <- steiger_filtering(dat3.3)
dat3.3 <- dat3.3[dat3.3$steiger_dir,]
res3.3 <- mr(dat3.3,method_list = c('mr_wald_ratio','mr_ivw'))
res3.3$exp <- 'Fenland'

res_ckd2 <- rbind(res1.3,res2.3,res3.3)

res_ckd2$out <- recode(res_ckd2$id.outcome,'ebi-a-GCST90018822'='UKBiobank+FinnGen','ebi-a-GCST003374'='CKDGen(trans-ancestry,2016)')

res_ckd <- rbind(res_ckd1,res_ckd2)

res_ckd1 <- res_ckd[res_ckd$out=='CKDGen(trans-ancestry,2019)',]
res_ckd2 <- res_ckd[res_ckd$out=='CKDGen(trans-ancestry,2016)',]
res_ckd3 <- res_ckd[res_ckd$out=='UKBiobank+FinnGen',]

res_ckd1.1 <- NULL
for(i in unique(res_ckd1$exposure)){
  dt <- res_ckd1[res_ckd1$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_ckd1.1 <- rbind(res_ckd1.1,dt)
}

res_ckd2.1 <- NULL
for(i in unique(res_ckd2$exposure)){
  dt <- res_ckd2[res_ckd2$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_ckd2.1 <- rbind(res_ckd2.1,dt)
}

res_ckd3.1 <- NULL
for(i in unique(res_ckd3$exposure)){
  dt <- res_ckd3[res_ckd3$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_ckd3.1 <- rbind(res_ckd3.1,dt)
}

###
res_ckd <- rbind(res_ckd1.1,res_ckd2.1,res_ckd3.1)
res_ckd <- generate_odds_ratios(res_ckd)
res_ckd$CI <- paste0(format(round(res_ckd$or,2),nsmall = 2),'(',format(round(res_ckd$or_lci95,2),nsmall = 2),'-',format(round(res_ckd$or_uci95,2),nsmall = 2),')')

###############################################gene ontology enrichment analysis##################################################################

gene1 <- c('C4A','FGF5','TCEA2','BTN3A2','GNPTG','IGFBP5','AGER','GCKR','YOD1')
gene2 <- c('GATM','SDCCAG8','BTN3A3','IDI2','HLA-DQA2','C4B','MFAP4','APOA4','DNAJC10','C2CD2L',
           'CEP170','PFKFB2','INHBA','INHBC','NFATC1','AIF1L','MICB','PLD3','AIF1','UMOD','HLA-E','LEAP2','GMPR')
map_dt <- bitr(c(gene1,gene2), fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
GO_database <- 'org.Hs.eg.db' 

GO <- enrichGO(map_dt$ENTREZID,
               OrgDb = GO_database,
               keyType = "ENTREZID",
               ont = "ALL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = T)
go_table <- as.data.frame(GO)

write.xlsx(go_table,"go_table.xlsx")
