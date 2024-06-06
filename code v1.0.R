
library(TwoSampleMR)
library(data.table)
library(openxlsx)
library(dplyr)
library(fdrtool)
library(meta)

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

#######################################e.g. read and format a pqtl data ######################################################

pqtl3 <- read.xlsx('pqtl3.xlsx')

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
  min_pval=0
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

res <- rbind(res1,res2,res3)
fdr <- fdrtool(res$pval,statistic="pvalue",plot=F)
res$qval <- fdr$qval 
res.fdr <- res[res$qval<0.05,]

##############combined effects by Meta analysis

res1.sig <- res1[res1$exposure%in%res.fdr$exposure,]
res1.sig$group <- 'cis'

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

##################sensitivity analysis example: results with both cis+trans pqtl

dat1$exposure <- gsub(' (cis)','',dat1$exposure,fixed = T)
dat1$exposure <- gsub(' (trans)','',dat1$exposure,fixed = T)
dat1a <- dat1[dat1$exposure%in%res1.sig$exposure,]
dat1a$id.exposure <- dat1a$exposure
res1a <- mr(dat1a,method_list = c('mr_wald_ratio','mr_ivw'))
res1a$group <- 'cis_trans'

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
egger1 <- mr_pleiotropy_test(dat1a_cis)
hetero1 <- mr_heterogeneity(dat1a_cis)

###############################################TWAS results####################################################
################################################################################################################

#mapping all coding genes
pqtl1_gene <- read.xlsx('pqtl1.xlsx')
pqtl1_gene <- unique(pqtl1_gene[pqtl1_gene$target%in%res.fdr$exposure,'gene'])

pqtl2_gene <- read.xlsx('pqtl2.xlsx')
pqtl2_gene <- unique(pqtl2_gene[pqtl2_gene$target%in%res.fdr$exposure,'gene'])

pqtl3_gene <- read.xlsx('pqtl3.xlsx')
pqtl3_gene <- unique(pqtl3_gene[pqtl3_gene$Target%in%res.fdr$exposure,'gene'])

all_gene <- c(pqtl1_gene,pqtl2_gene,pqtl3_gene)
all_gene <- unique(unlist(strsplit(all_gene,'|',fixed = T)))

gene_list <- read.table("2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",header = T)
gene_list <- unique(gene_list[,c('Gene','GeneSymbol')]) # acquire the GeneSymbol

################### eQTLGen results

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

#########################################example of sensitivity and replication analysis#######################

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
res1.2$out <- 'CKDGen(trans-ancestry,2019)'

res_ckd1 <- rbind(res1.2,res2.2,res3.2)

res_ckd1.1 <- NULL   #meta results
for(i in unique(res_ckd1$exposure)){
  dt <- res_ckd1[res_ckd1$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_ckd1.1 <- rbind(res_ckd1.1,dt)
}

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
