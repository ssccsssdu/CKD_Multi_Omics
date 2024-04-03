
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
  effect_allele_col ="Amin",
  other_allele_col = "Amaj",
  pval_col = "p",
  eaf_col = "maf",
  phenotype_col='exposure',
  min_pval=0,
  units_col='cis',
  chr_col = "Chr",
  pos_col = "Pos_b37"
)

pqtl1$type <- pqtl1$units.exposure
pqtl1$units.exposure <- 'SD'
pqtl1$samplesize.exposure <- 35597
pqtl1 <- pqtl1[-grep(',',pqtl1$SNP),]  #snp without snpid
pqtl1 <- pqtl1[-grep('chr',pqtl1$SNP),] #snp without snpid
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
  min_pval=0
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
  snp_col = "rsID",
  beta_col = "Effect",
  se_col = "SE",
  effect_allele_col ="EA",
  other_allele_col = "NEA",
  pval_col = "P.value",
  eaf_col = "EAF",
  phenotype_col='Target',
  units_col='cis.trans',
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

##########待定#########

dat1$exposure <- gsub(' (cis)','',dat1$exposure,fixed = T)
dat2$exposure <- gsub(' (cis)','',dat2$exposure,fixed = T)
dat3$exposure <- gsub(' (cis)','',dat3$exposure,fixed = T)

dat1.fdr <- dat1[dat1$exposure%in%res.fdr$exposure,]
dat2.fdr <- dat2[dat2$exposure%in%res.fdr$exposure,]
dat3.fdr <- dat3[dat3$exposure%in%res.fdr$exposure,]

dat1.fdr$dat <- 'Iceland'
dat2.fdr$dat <- 'UK Biobank'
dat3.fdr$dat <- 'Fenland'

dat.fdr <- dplyr::bind_rows(dat1.fdr,dat2.fdr,dat3.fdr)

#########Venn plot for overlap of proteins acorss three datasets

gg <- ggvenn(data=list(Iceland=res1$exposure,UKB=res2$exposure,Fenland=res3$exposure))

#########Volcano plot of proteinS associated with CKD

prot_all <- res[res$pval<0.05,]
prot_all <- prot_all%>%group_by(exposure)%>%mutate(count=n())
prot_all3 <- prot_all[prot_all$count==3,]
prot_all2 <- prot_all[prot_all$count==2,]
prot_all1 <- prot_all[prot_all$count==1,]

res$group <- ifelse(res$qval<0.05,'q<0.05(FDR corrected)','q>0.05(Insignificant)')
res$group <- factor(res$group,levels=c('q<0.05(FDR corrected)',"q>0.05(Insignificant)"))
res$label <- ifelse(res$qval<0.05,as.character(res$exposure),NA)
res$shape <- ifelse(res$exposure%in%prot_all3$exposure,'Repleted in 3 datasets',
                    ifelse(res$exposure%in%prot_all2$exposure,"Repleted in 2 datasets","Significant in 1 dataset"))
res$shape <- factor(res$shape,levels=c("Significant in 1 dataset",'Repleted in 2 datasets','Repleted in 3 datasets'))

res_Iceland <- res[res$dat=='Iceland',]

gg_pqtl1 <- ggplot(res_Iceland, aes(b, -log10(pval),shape=shape))+ 
  scale_x_continuous(limits=c(-0.9,0.6))+
  scale_y_continuous(limits=c(0,20))+
  geom_point(aes(color = group),alpha=0.6,size=4)+
  scale_color_manual(values = c('#ED0000FF','grey'))+
  geom_vline(xintercept = c(0), lty=3,color = 'black', lwd=0.5)+ 
  geom_hline(yintercept = -log10(5e-4), lty=3,color = 'darkred', lwd=1)+
  theme_bw()+
  theme(legend.title = element_blank(),legend.position=c(0.2,0.85),panel.grid = element_blank())+
  theme(legend.text = element_text(size=12),
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )+
  labs(title="Iceland",x = 'beta',y = '-log10 pvalue')+
  geom_text_repel(aes(x = b,y = -log10(pval),label=label),                       
                  max.overlaps = 100000,size=5,segment.color='black',show.legend=F)

res_UKB <- res[res$dat=='UK Biobank',]

gg_pqtl2 <- ggplot(res_UKB, aes(b, -log10(pval),shape=shape))+ 
  scale_x_continuous(limits=c(-0.7,1.1))+
  scale_y_continuous(limits=c(0,55))+
  geom_point(aes(color = group),alpha=0.6, size=4)+
  scale_color_manual(values = c('#ED0000FF','grey'))+
  geom_vline(xintercept = c(0), lty=3,color = 'black', lwd=0.5)+ 
  geom_hline(yintercept = -log10(5e-4), lty=3,color = 'darkred', lwd=1)+
  theme_bw()+
  theme(legend.title = element_blank(),legend.position="none",panel.grid = element_blank())+
  theme(legend.text = element_text(size=12),
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )+
  labs(title="UK Biobank",x = 'beta',y = '-log10 pvalue')+
  geom_text_repel(aes(x = b,y = -log10(pval),label=label),                       
                  max.overlaps = 100000,size=5,segment.color='black',show.legend=F)

res_Fenland <- res[res$dat=='Fenland',]

gg_pqtl3 <- ggplot(res_Fenland, aes(b, -log10(pval),shape=shape))+ 
  scale_x_continuous(limits=c(-0.7,0.6))+
  scale_y_continuous(limits=c(0,15))+
  geom_point(aes(color = group),alpha=0.6, size=4)+
  scale_color_manual(values = c('#ED0000FF','grey'))+
  geom_vline(xintercept = c(0), lty=3,color = 'black', lwd=0.5)+
  geom_hline(yintercept = -log10(5e-4), lty=3,color = 'darkred', lwd=1)+
  theme_bw()+
  theme(legend.title = element_blank(),legend.position="none",panel.grid = element_blank())+
  theme(legend.text = element_text(size=12),
    plot.title = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )+
  labs(title="Fenland",x = 'beta',y = '-log10 pvalue')+
  geom_text_repel(aes(x = b,y = -log10(pval),label=label),                       
                  max.overlaps = 100000,size=5,segment.color='black',show.legend=F)

CairoPDF('pqtl.pdf',width = 20, height = 8)
ggarrange(gg_pqtl1,gg_pqtl2,gg_pqtl3,nrow = 1,ncol = 3)
dev.off()

##############acquire combined effects by Meta analysis

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
res.sig$group <- factor(res.sig$data,levels = c('Iceland','UK Biobank','Fenland','Combined'))

##########forest plot of identified proteins############################################

##########defined the order according to beta

res.order <- NULL

for(i in unique(res.sig$exposure)){
  dt <- res.sig[res.sig$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt$source <- paste(dt$data,collapse = ',')
  dt <- dt[1,]
  dt$data <- 'Combined'
  res.order <- rbind(res.order,dt)
}

res.order <- res.order[order(res.order$b),]
level <- res.order$exposure

##########format res for forest plot

res.sig$exposure <- factor(res.sig$exposure,levels = level)
res.sig <- res.sig[order(res.sig$exposure,res.sig$data),]

res.sig <- generate_odds_ratios(res.sig)
res.sig$CI <- paste0(format(round(res.sig$or,2),nsmall = 2),'(',format(round(res.sig$or_lci95,2),nsmall = 2),'-',format(round(res.sig$or_uci95,2),nsmall = 2),')')

#write.xlsx(res.sig.adj,'res.sig.xlsx')

res.sig$CI <- ifelse(res.sig$num>1&res.sig$data!='Combined',NA,res.sig$CI)
res.sig$source <- ifelse(res.sig$num>1&res.sig$data!='Combined',NA,res.sig$data)

gg <- ggforestplot::forestplot(
  df = res.sig,
  name=exposure,
  estimate = b,
  se=se,
  pvalue = pval,
  psignif = 0.05,
  xlab = "OR(95%CI) of 1-SD increment in proteins on CKD",
  colour = data,
  logodds = T,
  xlim = c(-1.5, 20),
  xtickbreaks = c(0.5,0.8,1,1.5,2.0))+
  scale_color_manual(values = c('darkred','#42B540FF','#00468BFF','#925E9FFF'))+
  #ggforce::facet_col(facets = ~data,scales = "free_y",space = "free")+
  theme(legend.title = element_blank(),legend.position="top")+
  geom_text(aes(x = 10,label=CI))+
  geom_text(aes(x = 15,label=source))

CairoPDF('forest.pdf',width = 14, height = 12)
gg
dev.off()

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
res.a$group <- factor(res.a$data,levels = c('Iceland','UK Biobank','Fenland','Combined'))

res.a$exposure <- factor(res.a$exposure,levels = level)
res.a <- res.a[order(res.a$exposure,res.a$data),]

res.a <- generate_odds_ratios(res.a)
res.a$CI <- paste0(format(round(res.a$or,2),nsmall = 2),'(',format(round(res.a$or_lci95,2),nsmall = 2),'-',format(round(res.a$or_uci95,2),nsmall = 2),')')
res.a$CI <- ifelse(res.a$num>1&res.a$data!='Combined',NA,res.a$CI)
res.a$source <- ifelse(res.a$num>1&res.a$data!='Combined',NA,res.a$data)

gg_cis_trans <- ggforestplot::forestplot(
  df = res.a,
  name=exposure,
  estimate = b,
  se=se,
  pvalue = pval,
  psignif = 0.05,
  xlab = "OR(95%CI) of 1-SD increment in proteins on CKD",
  colour = data,
  logodds = T,
  xlim = c(-1.5, 20),
  xtickbreaks = c(0.5,0.8,1,1.5,2.0))+
  scale_color_manual(values = c('darkred','#42B540FF','#00468BFF','#925E9FFF'))+
  #ggforce::facet_col(facets = ~data,scales = "free_y",space = "free")+
  theme(legend.title = element_blank(),legend.position="top")+
  geom_text(aes(x = 10,label=CI))+
  geom_text(aes(x = 15,label=source))

CairoPDF('forest.cis_trans.pdf',width = 14, height = 12)
gg_cis_trans
dev.off()

########################################F statistics#####################################################

dat <- dat1a[dat1a$type=='cis',]
exp_out <- unique(dat[,c('exposure','outcome')])
R.square1 <- NULL
for(i in 1:nrow(exp_out)){
  SNP_all <- dat[dat$exposure==exp_out[i,1]&dat$outcome==exp_out[i,2],]
  Fvalue <- sum(SNP_all$beta.exposure^2/SNP_all$se.exposure^2,na.rm = T)
  Fvalue <- ifelse(Fvalue==0,NA,Fvalue)
  nsnp <- length(unique(SNP_all$SNP))
  N <- max(SNP_all$samplesize.exposure)
  R.s <- data.frame(exposure=exp_out[i,1],outcome=exp_out[i,2],Fvalue=Fvalue,SNP.keep=nsnp,n=N)
  R.square1 <- rbind(R.square1,R.s)
}

R.square1$dat <- 'Iceland'

dat <- dat2a[dat2a$type=='cis',]
exp_out <- unique(dat[,c('exposure','outcome')])
R.square2 <- NULL
for(i in 1:nrow(exp_out)){
  SNP_all <- dat[dat$exposure==exp_out[i,1]&dat$outcome==exp_out[i,2],]
  Fvalue <- sum(SNP_all$beta.exposure^2/SNP_all$se.exposure^2,na.rm = T)
  Fvalue <- ifelse(Fvalue==0,NA,Fvalue)
  nsnp <- length(unique(SNP_all$SNP))
  N <- max(SNP_all$samplesize.exposure)
  R.s <- data.frame(exposure=exp_out[i,1],outcome=exp_out[i,2],Fvalue=Fvalue,SNP.keep=nsnp,n=N)
  R.square2 <- rbind(R.square2,R.s)
}

R.square2$dat <- 'UK Biobank'

dat <- dat3a[dat3a$type=='cis',]
exp_out <- unique(dat[,c('exposure','outcome')])
R.square3 <- NULL
for(i in 1:nrow(exp_out)){
  SNP_all <- dat[dat$exposure==exp_out[i,1]&dat$outcome==exp_out[i,2],]
  Fvalue <- sum(SNP_all$beta.exposure^2/SNP_all$se.exposure^2,na.rm = T)
  Fvalue <- ifelse(Fvalue==0,NA,Fvalue)
  nsnp <- length(unique(SNP_all$SNP))
  N <- max(SNP_all$samplesize.exposure)
  R.s <- data.frame(exposure=exp_out[i,1],outcome=exp_out[i,2],Fvalue=Fvalue,SNP.keep=nsnp,n=N)
  R.square3 <- rbind(R.square3,R.s)
}

R.square3$dat <- 'Fenland'
R.square <- rbind(R.square1,R.square2,R.square3)

write.xlsx(R.square,'R.square.xlsx')

########

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
gene_list <- unique(gene_list[,c('Gene','GeneSymbol')])

###################eQTLGen results

gene1 <- read.table("ckd.smr",header=T)
gene1 <- merge(gene1,gene_list,by='Gene')
gene1 <- gene1[gene1$GeneSymbol%in%all_gene,]
gene1$Gene <- gene1$GeneSymbol

###################forest plot for egtlgen results

gene1.1 <- gene1
gene1.1$b <- gene1.1$b_SMR
gene1.1$se <- gene1.1$se_SMR
gene1.1 <- generate_odds_ratios(gene1.1)
gene1.1$CI <- paste0(format(round(gene1.1$or,2),nsmall = 2),'(',format(round(gene1.1$or_lci95,2),nsmall = 2),'-',format(round(gene1.1$or_uci95,2),nsmall = 2),')')
gene1.1$group <- ifelse(gene1.1$p_HEIDI>0.01,'HEIDI test passed','HEIDI test failed')
gene1.1 <- gene1.1[order(gene1.1$b_SMR),]

gg <- ggforestplot::forestplot(
  df = gene1.1,
  name=GeneSymbol,
  estimate = b_SMR,
  se=se_SMR,
  pvalue = p_SMR,
  psignif = 0.05,
  xlab = "OR(95%CI) of 1-SD increment in gene expression on CKD",
  colour = group,
  logodds = T,
  xlim = c(-0.2, 4),
  xtickbreaks = c(0.7, 1, 1.5, 2))+
  scale_color_manual(values = c('#B2182B','#2166AC'))+
  theme(legend.title = element_blank(),legend.position="top")+
  geom_text(aes(x = 3,label=CI))

CairoPDF(paste0("forest.eqtl.pdf"),width = 10, height = 10)
gg
dev.off()

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
gene12345$b_SMR <- ifelse(gene12345$p_HEIDI>0.01,gene12345$b_SMR,NA) #not show results did not passed HEIDI test
gene12345$p_SMR <- ifelse(gene12345$p_HEIDI>0.01,gene12345$p_SMR,NA)

corr <- gene12345
corr$tissue <- factor(corr$tissue,levels=unique(corr$tissue))

order <- c('SDCCAG8','CEP170','AGER','C4B','AIF1L',
           'DNAJC10','YOD1','HLA-E','C2CD2L','MICB','GCKR','INHBA','APOA4','LEAP2','AIF1','PLD3','TCEA2',
           'GATM','GMPR','PFKFB2','IGFBP5','FGF5','IDI2','GNPTG','NFATC1',
           'C4A','BTN3A3','BTN3A2','HLA-DQA2')

cor.mat <- dcast(corr,Gene~tissue,value.var='b_SMR')
row.names(cor.mat) <- cor.mat[,1]
cor.mat <- cor.mat[,-1]
cor.mat <- cor.mat[order,]

p.mat <- dcast(corr,Gene~tissue,value.var='p_SMR')
row.names(p.mat) <- p.mat[,1]
p.mat <- p.mat[,-1]
p.mat <- p.mat[order,]

cor.mat <- as.matrix(cor.mat)
p.mat <- as.matrix(p.mat)

cor.mat <- cbind(cor.mat[,c(1:3,53,4:52)])
p.mat <- cbind(p.mat[,c(1:3,53,4:52)])


colname <- colnames(cor.mat)
#write.csv(colname,"colname.csv",row.names = F)
## add sample size to colname
colname <- fread("colname.csv")
colname$colname <- paste0(colname$x,',',colname$n)
colnames(cor.mat) <- colname$colname
colnames(p.mat) <- colname$colname

col <- colorRampPalette(colors =c("darkgreen","white","red"),space="Lab")


CairoPDF(paste0("heatmap.pdf"),width = 12, height = 8)

corrplot::corrplot(cor.mat,is.corr = FALSE,col = col(200),method = 'color', type = 'full', tl.col = "black",
                   tl.srt = 45,tl.cex=0.7, diag = T, p.mat = p.mat, rect.col='black',rect.lwd=0.1,
                   na.label.col = "grey",na.label='-',addgrid.col='black',cl.ratio=0.1,
                   sig.level = c(0.05/29, 0.05), insig = 'label_sig',outline=T,cl.cex=0.7,
                   pch.cex = 1)

dev.off()

#########################################sensitivity and replication analysis#######################
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
# res_ckd2$out <- ifelse(res_ckd2$id.outcome=='ebi-a-GCST90018822','UKBiobank+FinnGen',NA)
# res_ckd2$out <- ifelse(res_ckd2$id.outcome=='ebi-a-GCST003374','CKDGen(trans-ancestry,2016)',res_ckd3$out)

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

res_ckd$exposure <- factor(res_ckd$exposure,levels = level)
res_ckd <- res_ckd[order(res_ckd$exposure),]

gg2 <- ggforestplot::forestplot(
  df = res_ckd,
  name=exposure,
  estimate = b,
  se=se,
  pvalue = pval,
  psignif = 0.05,
  xlab = "OR(95%CI) of 1-SD increment in proteins on CKD",
  colour = out,
  logodds = T,
  xlim = c(-1.5, 4),
  xtickbreaks = c(0.5, 1,1.5,2.0))+
  theme(legend.title = element_blank(),legend.position="top")

#########################renal function

ckdi25 <- fread("CKDGen_ckdi25_overall.txt")
ckdi25 <- ckdi25[ckdi25$RSID%in%c(pqtl1$SNP,pqtl2$SNP,pqtl3$SNP)]
ckdi25$p <- ckdi25$`P-value`
ckdi25$b <- log(ckdi25$OR)

ckdi25 <- format_data(
  ckdi25,
  type='outcome',
  snp_col = "RSID",
  beta_col = "b",
  se_col = "StdErr",
  effect_allele_col ="Allele1",
  other_allele_col = "Allele2",
  pval_col = "p",
  eaf_col = "Freq1",
  samplesize_col = "n_total_sum",
  units_col='cis',
  min_pval=0
)

ckdi25$outcome <- 'CKDi25'
ckdi25$id.outcome <- 'CKDi25'
ckdi25$units.outcome <- 'log odds'
ckdi25$ncase.outcome <- 19901
ckdi25$ncontrol.outcome <- 175244
ckdi25$samplesize.outcome <- ckdi25$ncase.outcome+ckdi25$ncontrol.outcome
ckdi25$prevalence.outcome <- ckdi25$ncase.outcome/ckdi25$samplesize.outcome

rapid3 <- fread("CKDGen_rapid3_overall.txt")
rapid3 <- rapid3[rapid3$RSID%in%c(pqtl1$SNP,pqtl2$SNP,pqtl3$SNP)]
rapid3$p <- rapid3$`P-value`
rapid3$b <- log(rapid3$OR)

rapid3 <- format_data(
  rapid3,
  type='outcome',
  snp_col = "RSID",
  beta_col = "b",
  se_col = "StdErr",
  effect_allele_col ="Allele1",
  other_allele_col = "Allele2",
  pval_col = "p",
  eaf_col = "Freq1",
  samplesize_col = "n_total_sum",
  units_col='cis',
  min_pval=0
)

rapid3$outcome <- 'Rapid3'
rapid3$id.outcome <- 'Rapid3'
rapid3$units.outcome <- 'log odds'
rapid3$ncase.outcome <-  34874
rapid3$ncontrol.outcome <- 107090
rapid3$samplesize.outcome <- rapid3$ncase.outcome+rapid3$ncontrol.outcome
rapid3$prevalence.outcome <- rapid3$ncase.outcome/rapid3$samplesize.outcome

eGFRcrea <- fread("metal_eGFR_meta1.TBL.map.annot.gc")
eGFRcrea <- eGFRcrea[eGFRcrea$RSID%in%c(pqtl1$SNP,pqtl2$SNP,pqtl3$SNP)]
eGFRcrea$p <- eGFRcrea$`P-value`

eGFRcrea <- format_data(
  eGFRcrea,
  type='outcome',
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col ="Allele1",
  other_allele_col = "Allele2",
  pval_col = "p",
  eaf_col = "Freq1",
  samplesize_col = "mac",
  units_col='cis',
  min_pval=0
)

eGFRcrea$outcome <- 'log(eGFRcrea)'
eGFRcrea$id.outcome <- 'log(eGFRcrea)'
eGFRcrea$units.outcome <- 'SD'
eGFRcrea$samplesize.outcome <- 1201909

eGFRcys <- fread("metal_eGFRcys_meta1.TBL.map.annot.gc")
eGFRcys <- eGFRcys[eGFRcys$RSID%in%c(pqtl1$SNP,pqtl2$SNP,pqtl3$SNP)]
eGFRcys$p <- eGFRcys$`P-value`

eGFRcys <- format_data(
  eGFRcys,
  type='outcome',
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col ="Allele1",
  other_allele_col = "Allele2",
  pval_col = "p",
  eaf_col = "Freq1",
  samplesize_col = "mac",
  units_col='cis',
  min_pval=0
)

eGFRcys$outcome <- 'log(eGFRcys)'
eGFRcys$id.outcome <- 'log(eGFRcys)'
eGFRcys$units.outcome <- 'SD'
eGFRcys$samplesize.outcome <- 1201909

kid_fun <- bind_rows(ckdi25,rapid3,eGFRcrea,eGFRcys)

#
dat1.4 <- harmonise_data(pqtl1a,kid_fun)
dat1.4 <- dat1.4[dat1.4$type=='cis',]
dat1.4$exposure <- gsub(' (cis)','',dat1.4$exposure,fixed = T)
dat1.4 <- dat1.4[dat1.4$exposure%in%res.fdr$exposure,]
dat1.4 <- dat1.4[dat1.4$mr_keep,]
dat1.4 <- dat1.4[!is.na(dat1.4$eaf.outcome),]
dat1.4 <- steiger_filtering(dat1.4)
dat1.4 <- dat1.4[dat1.4$steiger_dir,]
res1.4 <- mr(dat1.4,method_list = c('mr_wald_ratio','mr_ivw'))
res1.4$exp <- 'Iceland'

dat2.4 <- harmonise_data(pqtl2a,kid_fun)
dat2.4 <- dat2.4[dat2.4$type=='cis',]
dat2.4$exposure <- gsub(' (cis)','',dat2.4$exposure,fixed = T)
dat2.4 <- dat2.4[dat2.4$exposure%in%res.fdr$exposure,]
dat2.4 <- dat2.4[dat2.4$mr_keep,]
dat2.4 <- dat2.4[!is.na(dat2.4$eaf.outcome),]
dat2.4 <- steiger_filtering(dat2.4)
dat2.4 <- dat2.4[dat2.4$steiger_dir,]
res2.4 <- mr(dat2.4,method_list = c('mr_wald_ratio','mr_ivw'))
res2.4$exp <- 'UK Biobank'

dat3.4 <- harmonise_data(pqtl3a,kid_fun)
dat3.4 <- dat3.4[dat3.4$type=='cis',]
dat3.4$exposure <- gsub(' (cis)','',dat3.4$exposure,fixed = T)
dat3.4 <- dat3.4[dat3.4$exposure%in%res.fdr$exposure,]
dat3.4 <- dat3.4[dat3.4$mr_keep,]
dat3.4 <- dat3.4[!is.na(dat3.4$eaf.outcome),]
dat3.4 <- steiger_filtering(dat3.4)
dat3.4 <- dat3.4[dat3.4$steiger_dir,]
res3.4 <- mr(dat3.4,method_list = c('mr_wald_ratio','mr_ivw'))
res3.4$exp <- 'Fenland'

res_kid_fun <- rbind(res1.4,res2.4,res3.4)

res_kid_fun1 <- res_kid_fun[res_kid_fun$outcome=='CKDi25',]
res_kid_fun2 <- res_kid_fun[res_kid_fun$outcome=='Rapid3',]
res_kid_fun3 <- res_kid_fun[res_kid_fun$outcome=='log(eGFRcrea)',]
res_kid_fun4 <- res_kid_fun[res_kid_fun$outcome=='log(eGFRcys)',]

res_fun1.1 <- NULL
for(i in unique(res_kid_fun1$exposure)){
  dt <- res_kid_fun1[res_kid_fun1$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_fun1.1 <- rbind(res_fun1.1,dt)
}

res_fun2.1 <- NULL
for(i in unique(res_kid_fun2$exposure)){
  dt <- res_kid_fun2[res_kid_fun2$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_fun2.1 <- rbind(res_fun2.1,dt)
}

res_fun3.1 <- NULL
for(i in unique(res_kid_fun3$exposure)){
  dt <- res_kid_fun3[res_kid_fun3$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_fun3.1 <- rbind(res_fun3.1,dt)
}

res_fun4.1 <- NULL
for(i in unique(res_kid_fun4$exposure)){
  dt <- res_kid_fun4[res_kid_fun4$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_fun4.1 <- rbind(res_fun4.1,dt)
}

###
res_kid_fun <- rbind(res_fun1.1,res_fun2.1,res_fun3.1,res_fun4.1)
res_kid_fun <- generate_odds_ratios(res_kid_fun)
res_kid_fun$CI <- paste0(format(round(res_kid_fun$or,2),nsmall = 2),'(',format(round(res_kid_fun$or_lci95,2),nsmall = 2),'-',format(round(res_kid_fun$or_uci95,2),nsmall = 2),')')

res_kid_fun1 <- res_kid_fun[res_kid_fun$outcome=='CKDi25',]
res_kid_fun2 <- res_kid_fun[res_kid_fun$outcome=='Rapid3',]
res_kid_fun3 <- res_kid_fun[res_kid_fun$outcome=='log(eGFRcrea)',]
res_kid_fun4 <- res_kid_fun[res_kid_fun$outcome=='log(eGFRcys)',]
res_kid_fun3$b <- -res_kid_fun3$b
res_kid_fun3$CI2 <- paste0(format(round(res_kid_fun3$b,3),nsmall = 3),'(',format(round(res_kid_fun3$b-1.96*res_kid_fun3$se,3),nsmall = 3),'~',format(round(res_kid_fun3$b+1.96*res_kid_fun3$se,3),nsmall = 3),')')
res_kid_fun4$b <- -res_kid_fun4$b
res_kid_fun4$CI2 <- paste0(format(round(res_kid_fun4$b,3),nsmall = 3),'(',format(round(res_kid_fun4$b-1.96*res_kid_fun4$se,3),nsmall = 3),'~',format(round(res_kid_fun4$b+1.96*res_kid_fun4$se,3),nsmall = 3),')')

res_kid_fun12 <- rbind(res_kid_fun1,res_kid_fun2)
res_kid_fun12$exposure <- factor(res_kid_fun12$exposure,levels = level)
res_kid_fun12 <- res_kid_fun12[order(res_kid_fun12$exposure),]

gg3.1 <- ggforestplot::forestplot(
  df = res_kid_fun12,
  name=exposure,
  estimate = b,
  se=se,
  pvalue = pval,
  psignif = 0.05,
  xlab = "OR(95%CI) of 1-SD increment in proteins on rapid kidney function decline",
  colour = outcome,
  logodds = F,
  xlim = c(-0.5, 1.2),
  xtickbreaks = c(0, 0.5,1,1.5))+
  theme(legend.title = element_blank(),legend.position="top")

res_kid_fun34 <- rbind(res_kid_fun3,res_kid_fun4)
res_kid_fun34$exposure <- factor(res_kid_fun34$exposure,levels = level)
res_kid_fun34 <- res_kid_fun34[order(res_kid_fun34$exposure),]

gg3.2 <- ggforestplot::forestplot(
  df = res_kid_fun34,
  name=exposure,
  estimate = b,
  se=se,
  pvalue = pval,
  psignif = 0.05,
  xlab = "β(95%CI) of 1-SD increment in proteins on decreased eGFR",
  colour = outcome,
  logodds = F,
  xlim = c(-0.08, 0.1),
  xtickbreaks = c(-0.05, 0.05))+
  theme(legend.title = element_blank(),legend.position="top")

res_ckd_table <- bind_rows(res_ckd,res_kid_fun12,res_kid_fun34)
#write.xlsx(res_ckd_table,"res_ckd_table.xlsx")

###########################CKD clinical types

pqtl1a <- pqtl1
pqtl1a <- pqtl1a[pqtl1a$type=='cis',]
pqtl1a$exposure <- gsub(' (cis)','',pqtl1a$exposure,fixed = T)
pqtl1a <- pqtl1a[pqtl1a$exposure%in%res.fdr$exposure,]
ckd4 <- extract_outcome_data(snps=pqtl1a$SNP,outcomes=c(
  'ebi-a-GCST90018866', #IgA nephropathy
  'ebi-a-GCST010005', #Membranous nephropathy
  'ebi-a-GCST90018820', #Chronic glomerulonephritis
  'ebi-a-GCST90018884',#Nephrotic syndrome
  "ebi-a-GCST90018832", #Diabetic nephropathy,
  "finn-b-N14_CHRONTUBULOINTNEPHRITIS"
))
dat1.4 <- harmonise_data(pqtl1a,ckd4)
dat1.4 <- dat1.4[dat1.4$mr_keep,]
dat1.4$eaf.outcome <- ifelse(is.na(dat1.4$eaf.outcome),dat1.4$eaf.exposure,dat1.4$eaf.outcome)
dat1.4 <- dat1.4[!is.na(dat1.4$eaf.outcome),]
dat1.4$units.outcome <- 'log odds'
dat1.4$ncase.outcome <- NA
dat1.4$ncontrol.outcome <- NA
dat1.4$ncase.outcome <- ifelse(dat1.4$id.outcome=='ebi-a-GCST90018866',15587,dat1.4$ncase.outcome)
dat1.4$ncontrol.outcome <- ifelse(dat1.4$id.outcome=='ebi-a-GCST90018866',462197,dat1.4$ncontrol.outcome)
dat1.4$ncase.outcome <- ifelse(dat1.4$id.outcome=='ebi-a-GCST010005',	2150,dat1.4$ncase.outcome)
dat1.4$ncontrol.outcome <- ifelse(dat1.4$id.outcome=='ebi-a-GCST010005',5829,dat1.4$ncontrol.outcome)
dat1.4$ncase.outcome <- ifelse(dat1.4$id.outcome=='ebi-a-GCST90018820',	566,dat1.4$ncase.outcome)
dat1.4$ncontrol.outcome <- ifelse(dat1.4$id.outcome=='ebi-a-GCST90018820',475255,dat1.4$ncontrol.outcome)
dat1.4$ncase.outcome <- ifelse(dat1.4$id.outcome=='ebi-a-GCST90018884',	775,dat1.4$ncase.outcome)
dat1.4$ncontrol.outcome <- ifelse(dat1.4$id.outcome=='ebi-a-GCST90018884',475255,dat1.4$ncontrol.outcome)
dat1.4$ncase.outcome <- ifelse(dat1.4$id.outcome=='ebi-a-GCST90018832',	1032,dat1.4$ncase.outcome)
dat1.4$ncontrol.outcome <- ifelse(dat1.4$id.outcome=='ebi-a-GCST90018832',451248,dat1.4$ncontrol.outcome)
dat1.4$ncase.outcome <- ifelse(dat1.4$id.outcome=='finn-b-N14_CHRONTUBULOINTNEPHRITIS',620,dat1.4$ncase.outcome)
dat1.4$ncontrol.outcome <- ifelse(dat1.4$id.outcome=='finn-b-N14_CHRONTUBULOINTNEPHRITIS',201028,dat1.4$ncontrol.outcome)

dat1.4$samplesize.outcome <- dat1.4$ncase.outcome+dat1.4$ncontrol.outcome
dat1.4$prevalence.outcome <- dat1.4$ncase.outcome/dat1.4$samplesize.outcome
dat1.4 <- steiger_filtering(dat1.4)
dat1.4 <- dat1.4[dat1.4$steiger_dir,]
res1.4 <- mr(dat1.4,method_list = c('mr_wald_ratio','mr_ivw'))
res1.4$exp <- 'Iceland'

#############################

pqtl2a <- pqtl2
pqtl2a <- pqtl2a[pqtl2a$type=='cis',]
pqtl2a$exposure <- gsub(' (cis)','',pqtl2a$exposure,fixed = T)
pqtl2a <- pqtl2a[pqtl2a$exposure%in%res.fdr$exposure,]
ckd4 <- extract_outcome_data(snps=pqtl2a$SNP,outcomes=c(
  'ebi-a-GCST90018866', #IgA nephropathy
  'ebi-a-GCST010005', #Membranous nephropathy
  'ebi-a-GCST90018820', #Chronic glomerulonephritis
  'ebi-a-GCST90018884',#Nephrotic syndrome
  "ebi-a-GCST90018832", #Diabetic nephropathy
  "finn-b-N14_CHRONTUBULOINTNEPHRITIS"
))

dat2.4 <- harmonise_data(pqtl2a,ckd4)
dat2.4 <- dat2.4[dat2.4$mr_keep,]
dat2.4$eaf.outcome <- ifelse(is.na(dat2.4$eaf.outcome),dat2.4$eaf.exposure,dat2.4$eaf.outcome)
dat2.4 <- dat2.4[!is.na(dat2.4$eaf.outcome),]
dat2.4$units.outcome <- 'log odds'
dat2.4$ncase.outcome <- NA
dat2.4$ncontrol.outcome <- NA
dat2.4$ncase.outcome <- ifelse(dat2.4$id.outcome=='ebi-a-GCST90018866',15587,dat2.4$ncase.outcome)
dat2.4$ncontrol.outcome <- ifelse(dat2.4$id.outcome=='ebi-a-GCST90018866',462197,dat2.4$ncontrol.outcome)
dat2.4$ncase.outcome <- ifelse(dat2.4$id.outcome=='ebi-a-GCST010005',	2150,dat2.4$ncase.outcome)
dat2.4$ncontrol.outcome <- ifelse(dat2.4$id.outcome=='ebi-a-GCST010005',5829,dat2.4$ncontrol.outcome)
dat2.4$ncase.outcome <- ifelse(dat2.4$id.outcome=='ebi-a-GCST90018820',	566,dat2.4$ncase.outcome)
dat2.4$ncontrol.outcome <- ifelse(dat2.4$id.outcome=='ebi-a-GCST90018820',475255,dat2.4$ncontrol.outcome)
dat2.4$ncase.outcome <- ifelse(dat2.4$id.outcome=='ebi-a-GCST90018884',	775,dat2.4$ncase.outcome)
dat2.4$ncontrol.outcome <- ifelse(dat2.4$id.outcome=='ebi-a-GCST90018884',475255,dat2.4$ncontrol.outcome)
dat2.4$ncase.outcome <- ifelse(dat2.4$id.outcome=='ebi-a-GCST90018832',	1032,dat2.4$ncase.outcome)
dat2.4$ncontrol.outcome <- ifelse(dat2.4$id.outcome=='ebi-a-GCST90018832',451248,dat2.4$ncontrol.outcome)
dat2.4$ncase.outcome <- ifelse(dat2.4$id.outcome=='finn-b-N14_CHRONTUBULOINTNEPHRITIS',620,dat2.4$ncase.outcome)
dat2.4$ncontrol.outcome <- ifelse(dat2.4$id.outcome=='finn-b-N14_CHRONTUBULOINTNEPHRITIS',201028,dat2.4$ncontrol.outcome)

dat2.4$samplesize.outcome <- dat2.4$ncase.outcome+dat2.4$ncontrol.outcome
dat2.4$prevalence.outcome <- dat2.4$ncase.outcome/dat2.4$samplesize.outcome
dat2.4 <- steiger_filtering(dat2.4)
dat2.4 <- dat2.4[dat2.4$steiger_dir,]
res2.4 <- mr(dat2.4,method_list = c('mr_wald_ratio','mr_ivw'))
res2.4$exp <- 'UK Biobank'

#####################

pqtl3a <- pqtl3
pqtl3a <- pqtl3a[pqtl3a$type=='cis',]
pqtl3a$exposure <- gsub(' (cis)','',pqtl3a$exposure,fixed = T)
pqtl3a <- pqtl3a[pqtl3a$exposure%in%res.fdr$exposure,]
ckd4 <- extract_outcome_data(snps=pqtl3a$SNP,outcomes=c(
  'ebi-a-GCST90018866', #IgA nephropathy
  'ebi-a-GCST010005', #Membranous nephropathy
  'ebi-a-GCST90018820', #Chronic glomerulonephritis
  'ebi-a-GCST90018884',#Nephrotic syndrome
  "ebi-a-GCST90018832", #Diabetic nephropathy
  "finn-b-N14_CHRONTUBULOINTNEPHRITIS"
))
dat3.4 <- harmonise_data(pqtl3a,ckd4)
dat3.4 <- dat3.4[dat3.4$mr_keep,]
dat3.4$eaf.outcome <- ifelse(is.na(dat3.4$eaf.outcome),dat3.4$eaf.exposure,dat3.4$eaf.outcome)
dat3.4 <- dat3.4[!is.na(dat3.4$eaf.outcome),]
dat3.4$units.outcome <- 'log odds'
dat3.4$ncase.outcome <- NA
dat3.4$ncontrol.outcome <- NA
dat3.4$ncase.outcome <- ifelse(dat3.4$id.outcome=='ebi-a-GCST90018866',15587,dat3.4$ncase.outcome)
dat3.4$ncontrol.outcome <- ifelse(dat3.4$id.outcome=='ebi-a-GCST90018866',462197,dat3.4$ncontrol.outcome)
dat3.4$ncase.outcome <- ifelse(dat3.4$id.outcome=='ebi-a-GCST010005',	2150,dat3.4$ncase.outcome)
dat3.4$ncontrol.outcome <- ifelse(dat3.4$id.outcome=='ebi-a-GCST010005',5829,dat3.4$ncontrol.outcome)
dat3.4$ncase.outcome <- ifelse(dat3.4$id.outcome=='ebi-a-GCST90018820',	566,dat3.4$ncase.outcome)
dat3.4$ncontrol.outcome <- ifelse(dat3.4$id.outcome=='ebi-a-GCST90018820',475255,dat3.4$ncontrol.outcome)
dat3.4$ncase.outcome <- ifelse(dat3.4$id.outcome=='ebi-a-GCST90018884',	775,dat3.4$ncase.outcome)
dat3.4$ncontrol.outcome <- ifelse(dat3.4$id.outcome=='ebi-a-GCST90018884',475255,dat3.4$ncontrol.outcome)
dat3.4$ncase.outcome <- ifelse(dat3.4$id.outcome=='ebi-a-GCST90018832',	1032,dat3.4$ncase.outcome)
dat3.4$ncontrol.outcome <- ifelse(dat3.4$id.outcome=='ebi-a-GCST90018832',451248,dat3.4$ncontrol.outcome)
dat3.4$ncase.outcome <- ifelse(dat3.4$id.outcome=='finn-b-N14_CHRONTUBULOINTNEPHRITIS',620,dat3.4$ncase.outcome)
dat3.4$ncontrol.outcome <- ifelse(dat3.4$id.outcome=='finn-b-N14_CHRONTUBULOINTNEPHRITIS',201028,dat3.4$ncontrol.outcome)

dat3.4$samplesize.outcome <- dat3.4$ncase.outcome+dat3.4$ncontrol.outcome
dat3.4$prevalence.outcome <- dat3.4$ncase.outcome/dat3.4$samplesize.outcome
dat3.4 <- steiger_filtering(dat3.4)
dat3.4 <- dat3.4[dat3.4$steiger_dir,]
res3.4 <- mr(dat3.4,method_list = c('mr_wald_ratio','mr_ivw'))
res3.4$exp <- 'Fenland'

res_ckd4 <- rbind(res1.4,res2.4,res3.4)
res_ckd4$out <- ifelse(res_ckd4$id.outcome=='finn-b-N14_CHRONTUBULOINTNEPHRITIS','Chronic tubulo-interstitial nephritis',NA)
res_ckd4$out <- ifelse(res_ckd4$id.outcome=='ebi-a-GCST90018866','IgA nephropathy',res_ckd4$out)
res_ckd4$out <- ifelse(res_ckd4$id.outcome=='ebi-a-GCST010005','Membranous nephropathy',res_ckd4$out)
res_ckd4$out <- ifelse(res_ckd4$id.outcome=='ebi-a-GCST90018820','Chronic glomerulonephritis',res_ckd4$out)
res_ckd4$out <- ifelse(res_ckd4$id.outcome=='ebi-a-GCST90018884','Nephrotic syndrome',res_ckd4$out)
res_ckd4$out <- ifelse(res_ckd4$id.outcome=='ebi-a-GCST90018832','Diabetic nephropathy',res_ckd4$out)

res_ckd4.1 <- res_ckd4[res_ckd4$out=='Chronic tubulo-interstitial nephritis',]
res_ckd4.2 <- res_ckd4[res_ckd4$out=='IgA nephropathy',]
res_ckd4.3 <- res_ckd4[res_ckd4$out=='Membranous nephropathy',]
res_ckd4.4 <- res_ckd4[res_ckd4$out=='Chronic glomerulonephritis',]
res_ckd4.5 <- res_ckd4[res_ckd4$out=='Nephrotic syndrome',]
res_ckd4.6 <- res_ckd4[res_ckd4$out=='Diabetic nephropathy',]

res_ckd4.11 <- NULL
for(i in unique(res_ckd4.1$exposure)){
  dt <- res_ckd4.1[res_ckd4.1$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_ckd4.11 <- rbind(res_ckd4.11,dt)
}

res_ckd4.21 <- NULL
for(i in unique(res_ckd4.2$exposure)){
  dt <- res_ckd4.2[res_ckd4.2$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_ckd4.21 <- rbind(res_ckd4.21,dt)
}

res_ckd4.31 <- NULL
for(i in unique(res_ckd4.3$exposure)){
  dt <- res_ckd4.3[res_ckd4.3$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_ckd4.31 <- rbind(res_ckd4.31,dt)
}

res_ckd4.41 <- NULL
for(i in unique(res_ckd4.4$exposure)){
  dt <- res_ckd4.4[res_ckd4.4$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_ckd4.41 <- rbind(res_ckd4.41,dt)
}

res_ckd4.51 <- NULL
for(i in unique(res_ckd4.5$exposure)){
  dt <- res_ckd4.5[res_ckd4.5$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_ckd4.51 <- rbind(res_ckd4.51,dt)
}

res_ckd4.61 <- NULL
for(i in unique(res_ckd4.6$exposure)){
  dt <- res_ckd4.6[res_ckd4.6$exposure==i,]
  met <- metagen(b,se,sm="OR",data=dt)
  dt$b <- met$TE.fixed
  dt$se <- met$seTE.fixed
  dt$pval <- met$pval.fixed
  dt <- dt[1,]
  res_ckd4.61 <- rbind(res_ckd4.61,dt)
}

###
res_ckd <- rbind(res_ckd4.11,res_ckd4.21,res_ckd4.31,res_ckd4.41,res_ckd4.51,res_ckd4.61)
res_ckd <- generate_odds_ratios(res_ckd)
res_ckd$CI <- paste0(format(round(res_ckd$or,2),nsmall = 2),'(',format(round(res_ckd$or_lci95,2),nsmall = 2),'-',format(round(res_ckd$or_uci95,2),nsmall = 2),')')

res_ckd$exposure <- factor(res_ckd$exposure,levels = level)
res_ckd <- res_ckd[order(res_ckd$exposure),]
#write.xlsx(res_ckd16,"res_ckd16.xlsx")

res_ckd.sig <- res_ckd[res_ckd$pval<0.05,]

gg4 <- ggforestplot::forestplot(
  df = res_ckd.sig,
  name=exposure,
  estimate = b,
  se=se,
  pvalue = pval,
  psignif = 0.05,
  xlab = "OR(95%CI) of 1-SD increment in proteins on six CKD clinical types",
  colour = out,
  logodds = T,
  xlim = c(-0.5, 10),
  xtickbreaks = c(0.5, 1,1.5))+
  theme(legend.title = element_blank(),legend.position="top")

CairoPDF(paste0("forest.pqtl4.pdf"),width = 18, height = 20)
ggarrange(gg2,gg3.2,gg3.1,gg4,nrow = 2,ncol = 2)
dev.off()

##################balloonplot of identified proteins on CKD-related phenotypes######################

dt1 <- res.order[,c('exposure','outcome','b','pval')]
dt1$outcome <- 'CKD(CKDGen,European)'
dt1 <- dt1[order(dt1$b),]

dt1.1 <- dt1[dt1$b<0,]
dt1.2 <- dt1[dt1$b>0,]

dt2 <- rbind(res_ckd1,res_ckd2,res_ckd3)
dt2 <- dt2[,c('exposure','out','b','pval')]
names(dt2)[2] <- 'outcome'

dt3 <- bind_rows(res_kid_fun1,res_kid_fun2,res_kid_fun3,res_kid_fun4)
dt3 <- dt3[,c('exposure','outcome','b','pval')]

dt4 <- bind_rows(res_ckd4.1,res_ckd4.2,res_ckd4.3,res_ckd4.4,res_ckd4.5,res_ckd4.6)
dt4 <- dt4[,c('exposure','out','b','pval')]
names(dt4)[2] <- 'outcome'

dt <- rbind(dt1,dt2,dt3,dt4)
dtt <- dt[dt$pval<0.05,]

dtt$b <- ifelse(dtt$outcome=='log(eGFRcrea)',10*dtt$b,dtt$b)
dtt$b <- ifelse(dtt$outcome=='log(eGFRcys)',10*dtt$b,dtt$b)

dtt$outcome <- recode(dtt$outcome,"UKBiobank+FinnGen"='CKD(UKBiobank+FinnGen)','CKDGen(trans-ancestry,2019)'='CKD(CKDGen trans-ancestry,2019)','CKDGen(trans-ancestry,2016)'='CKD(CKDGen trans-ancestry,2016)')
dtt$outcome <- ifelse(dtt$outcome=='log(eGFRcrea)','10-unit decreased log(eGFRcrea)',dtt$outcome)
dtt$outcome <- ifelse(dtt$outcome=='log(eGFRcys)','10-unit decreased log(eGFRcys)',dtt$outcome)
dtt$outcome <- ifelse(dtt$outcome=='ckdi25','CKDi25',dtt$outcome)
dtt$outcome <- ifelse(dtt$outcome=='rapid3','Rapid3',dtt$outcome)
dtt$outcome <- ifelse(dtt$outcome=='rapid3','Rapid3',dtt$outcome)

dtt$outcome <- factor(dtt$outcome,levels = unique(dtt$outcome))
dtt$exposure <- factor(dtt$exposure,levels = rev(dt1$exposure))

dtt$log10P <- -log10(dtt$pval)
dtt$group <- ifelse(dtt$exposure%in%dt1.1$exposure,'Favour to decreased risk','Favour to increased risk')

plott <- ggballoonplot(dtt,x='outcome',y='exposure',
                       size='log10P',
                       fill='b',
                       rotate.x.text=FALSE,
                       #facet.by='group',
                       ggtheme=scale_fill_gradient2(low = "green",mid='white', high = "red"))+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45,face = "bold",vjust=1,hjust=1))+
  theme(axis.text.y=element_text(face = "bold"))+
  facet_grid(group~.,scales = 'free_y')

CairoPDF(paste0("ggballoonplot.pdf"),width = 10, height = 8)
print(plott)
dev.off()

###############################################gene ontology enrichment analysis##################################################################

library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

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

GO2 <- pairwise_termsim(GO)
go2 <- enrichplot::emapplot(GO2,showCategory = 20, color = "p.adjust", layout = "kk",cex_line=0.3,cex_label_category=0.5)#通路间关联网络图

CairoPDF(paste0("Go.pdf"),width = 8, height = 8)
go2
dev.off()


######################################Colocalization analysis#########################################################

pqtl1_gene <- read.xlsx('deCODE_pqtl.xlsx')
pqtl2_gene <- read.xlsx('UKB_PPP_pqtl.xlsx')
pqtl3_gene <- read.xlsx('Fenland_pqtl.xlsx')

ckd <- fread("CKD_overall_EA_JW_20180223_nstud23.dbgap.txt")
ckd$p <- ckd$`P-value`
ckd <- ckd[ckd$Freq1>0&ckd$Freq1<1,]
ckd$Freq1 <- ifelse(ckd$Freq1>0.5,1-ckd$Freq1,ckd$Freq1)

file1 <- list.files("data_ukb",full.names = T)
file2 <- list.files("data_ukb")
file2 <- gsub('.txt','',file2)

coloc_res_ukb <- NULL

CairoPDF(paste0("coloc_ukb.pdf"),width = 6, height = 6)

for(i in 1:length(file2)){
  
  dat <- fread(file1[i])
  cis_SNP <- pqtl2_gene[pqtl2_gene$target==file2[i]&pqtl2_gene$cis=='cis',]$SNP
  chr <- unique(dat[dat$rsids%in%cis_SNP,]$Chrom)
  
  for(j in 1:length(cis_SNP)){
    
    pos1 <- min(dat[dat$rsids%in%cis_SNP[j],]$Pos)-500000
    pos2 <- max(dat[dat$rsids%in%cis_SNP[j],]$Pos)+500000
    dats <- dat[dat$Chrom==chr&between(dat$Pos,pos1,pos2),]
    dats <- as.data.frame(dats)
    dats <- dats[!is.na(dats$rsids),]
    dats <- dats[!is.na(dats$Beta),]
    dats <- dats[!duplicated(dats$rsids),]
    dats <- dats[dats$ImpFreqA1>0&dats$ImpFreqA1<1,]
    dats$ImpFreqA1 <- ifelse(dats$ImpFreqA1>0.5,1-dats$ImpFreqA1,dats$ImpFreqA1)
    ckd2 <- ckd[ckd$RSID%in%dats$rsids,]
    x1 <- list(beta=dats$Beta,varbeta=dats$SE^2,pvalues=dats$Pval,MAF=dats$ImpFreqA1,snp=dats$rsids,type='quant',N=54219)
    x2 <- list(beta=ckd2$Effect,varbeta=ckd2$StdErr^2,pvalues=ckd2$p,MAF=ckd2$Freq1,snp=ckd2$RSID,type='cc')
    res <- coloc.abf(x1, x2)
    res <- data.frame(nsnps=unlist(res$summary)[1],
                      PP.H0.abf=unlist(res$summary)[2],
                      PP.H1.abf=unlist(res$summary)[3],
                      PP.H2.abf=unlist(res$summary)[4],
                      PP.H3.abf=unlist(res$summary)[5],
                      PP.H4.abf=unlist(res$summary)[6],
                      protein=file2[i],
                      data='UKB-PPP',
                      snp=cis_SNP[j])
    
    coloc_res_ukb <- rbind(coloc_res_ukb,res)
    
    coloc_plot <- locusplotr::gg_locusplot(df=ckd2,
                                           lead_snp=cis_SNP[j],
                                           rsid=RSID,
                                           chrom=Chr,
                                           pos=Pos_b37,
                                           ref = Allele1,
                                           alt = Allele2,
                                           p_value = p,
                                           plot_genes = T,
                                           plot_title=paste0('UKB-PPP:',file2[i],',PP.H4=',round(res$PP.H4.abf,2)))
    
    print(coloc_plot)
    
  }
  
  cat(i,'\n')
}


dev.off()

write.xlsx(coloc_res_ukb,"D:\\01-Research\\CKD_drug\\res2\\coloc_res_ukb.xlsx")

###########Iceland

file1 <- list.files("data_ice",full.names = T)
file2 <- list.files("data_ice")
file2 <- gsub('.txt','',file2)

coloc_res_ice <- NULL

CairoPDF(paste0("coloc_ice.pdf"),width = 6, height = 6)

for(i in 1:length(file2)){
  
  dat <- fread(file1[i])
  cis_SNP <- pqtl1_gene[pqtl1_gene$target==file2[i]&pqtl1_gene$cis=='cis',]$SNP
  cis_SNP <- cis_SNP[grep('rs',cis_SNP)]
  chr <- unique(dat[dat$rsids%in%cis_SNP,]$Chrom)
  
  for(j in 1:length(cis_SNP)){
    
    pos1 <- min(dat[dat$rsids%in%cis_SNP[j],]$Pos)-500000
    pos2 <- max(dat[dat$rsids%in%cis_SNP[j],]$Pos)+500000
    dats <- dat[dat$Chrom==chr&between(dat$Pos,pos1,pos2),]
    dats <- as.data.frame(dats)
    dats <- dats[!is.na(dats$rsids),]
    dats <- dats[!is.na(dats$Beta),]
    dats <- dats[!duplicated(dats$rsids),]
    dats <- dats[dats$ImpMAF>0&dats$ImpMAF<1,]
    dats$ImpMAF <- ifelse(dats$ImpMAF>0.5,1-dats$ImpMAF,dats$ImpMAF)
    ckd2 <- ckd[ckd$RSID%in%dats$rsids,]
    x1 <- list(beta=dats$Beta,varbeta=dats$SE^2,pvalues=dats$Pval,MAF=dats$ImpMAF,snp=dats$rsids,type='quant',N=35597)
    x2 <- list(beta=ckd2$Effect,varbeta=ckd2$StdErr^2,pvalues=ckd2$p,MAF=ckd2$Freq1,snp=ckd2$RSID,type='cc')
    res <- coloc.abf(x1, x2)
    res <- data.frame(nsnps=unlist(res$summary)[1],
                      PP.H0.abf=unlist(res$summary)[2],
                      PP.H1.abf=unlist(res$summary)[3],
                      PP.H2.abf=unlist(res$summary)[4],
                      PP.H3.abf=unlist(res$summary)[5],
                      PP.H4.abf=unlist(res$summary)[6],
                      protein=file2[i],
                      data='Iceland',
                      snp=cis_SNP[j])
    
    coloc_res_ice <- rbind(coloc_res_ice,res)
    
    temp <- try(locusplotr::gg_locusplot(df=ckd2,
                                         lead_snp=cis_SNP[j],
                                         rsid=RSID,
                                         chrom=Chr,
                                         pos=Pos_b37,
                                         ref = Allele1,
                                         alt = Allele2,
                                         p_value = p,
                                         plot_genes = T,
                                         plot_title=paste0('Iceland:',file2[i],',PP.H4=',round(res$PP.H4.abf,2))),silent=FALSE)
    
    if(!'try-error' %in% class(temp)){
      coloc_plot <- locusplotr::gg_locusplot(df=ckd2,
                                             lead_snp=cis_SNP[j],
                                             rsid=RSID,
                                             chrom=Chr,
                                             pos=Pos_b37,
                                             ref = Allele1,
                                             alt = Allele2,
                                             p_value = p,
                                             plot_genes = T,
                                             plot_title=paste0('Iceland:',file2[i],',PP.H4=',round(res$PP.H4.abf,2)))
      
      print(coloc_plot)
    }
  }
  
  cat(i,'\n')
}


dev.off()

write.xlsx(coloc_res_ice,"coloc_res_ice.xlsx")

############################Fen

file1 <- list.files("data_ice",full.names = T)
file2 <- list.files("data_ice")
file2 <- gsub('.txt','',file2)

coloc_res_fen <- NULL

CairoPDF(paste0("coloc_fen.pdf"),width = 6, height = 6)

for(i in 1:length(file2)){
  
  dat <- fread(file1[i])
  cis_SNP <- pqtl3_gene[pqtl3_gene$Target==file2[i]&pqtl3_gene$cis=='cis',]$rsID
  chr <- unique(dat[dat$rsids%in%cis_SNP,]$Chrom)
  
  for(j in 1:length(cis_SNP)){
    
    pos1 <- min(dat[dat$rsids%in%cis_SNP[j],]$Pos)-500000
    pos2 <- max(dat[dat$rsids%in%cis_SNP[j],]$Pos)+500000
    dats <- dat[dat$Chrom==chr&between(dat$Pos,pos1,pos2),]
    dats <- as.data.frame(dats)
    dats <- dats[!is.na(dats$rsids),]
    dats <- dats[!is.na(dats$Beta),]
    dats <- dats[!duplicated(dats$rsids),]
    dats <- dats[dats$ImpFreqA1>0&dats$ImpFreqA1<1,]
    dats$ImpFreqA1 <- ifelse(dats$ImpFreqA1>0.5,1-dats$ImpFreqA1,dats$ImpFreqA1)
    ckd2 <- ckd[ckd$RSID%in%dats$rsids,]
    x1 <- list(beta=dats$Beta,varbeta=dats$SE^2,pvalues=dats$Pval,MAF=dats$ImpFreqA1,snp=dats$rsids,type='quant',N=54219)
    x2 <- list(beta=ckd2$Effect,varbeta=ckd2$StdErr^2,pvalues=ckd2$p,MAF=ckd2$Freq1,snp=ckd2$RSID,type='cc')
    res <- coloc.abf(x1, x2)
    res <- data.frame(nsnps=unlist(res$summary)[1],
                      PP.H0.abf=unlist(res$summary)[2],
                      PP.H1.abf=unlist(res$summary)[3],
                      PP.H2.abf=unlist(res$summary)[4],
                      PP.H3.abf=unlist(res$summary)[5],
                      PP.H4.abf=unlist(res$summary)[6],
                      protein=file2[i],
                      data='UKB-PPP',
                      snp=cis_SNP[j])
    
    coloc_res_fen <- rbind(coloc_res_fen,res)
    
    coloc_plot <- locusplotr::gg_locusplot(df=ckd2,
                                           lead_snp=cis_SNP[j],
                                           rsid=RSID,
                                           chrom=Chr,
                                           pos=Pos_b37,
                                           ref = Allele1,
                                           alt = Allele2,
                                           p_value = p,
                                           plot_genes = T,
                                           plot_title=paste0('Iceland:',file2[i],',PP.H4=',round(res$PP.H4.abf,2)))
    
    print(coloc_plot)
    
  }
  
  cat(i,'\n')
}

dev.off()

write.xlsx(coloc_res_fen,"coloc_res_fen.xlsx")


