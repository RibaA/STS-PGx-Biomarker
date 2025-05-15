#####################################################################
## packages
#####################################################################

library(qs)
library(PharmacoGx)
library(data.table)
library(magrittr)
library(GSVA)
library(fgsea)
library(reshape2)
#library(enrichR)

###################################################
## Read KEGG and HALLMARK pathways from msigdb
###################################################

hallmark_pathway <- gmtPathways("/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Data/h.all.v2024.1.Hs.symbols.gmt")
kegg_pathway <- gmtPathways("/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Data/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt")

dat <- list(KEGG_pathway = kegg_pathway, HALLMARK_pathway = hallmark_pathway)
qsave(dat, file="/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Data/kegg_hallmark_pathway.qs")

###################################################
## GSEA: STS gene association result and HALLMARK
###################################################

dat <- qread(file="/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/drug/gene_drug_assoc_sts_late_integrated_celligner.qs")
dat <- dat[!is.na(dat$Coef), ]
drug <- unique(dat$drug)

gsea_hallmark <- lapply(1:length(drug), function(i){
  
  df <- dat[dat$drug == drug[i], ] 
  #ranks <- df$Coef
  ranks <- sign(df$Coef) * (-log10(df$Pval))
  names(ranks) <-df$Gene
  ranks <- sort(ranks, decreasing = T)
  
  fgseaRes <- fgsea(hallmark_pathway, 
                    ranks, 
                    minSize=15, 
                    maxSize = 500)
  fgseaRes <- fgseaRes[!is.na(fgseaRes$padj), ]
  
  data.frame(drug = drug[i], 
             fgseaRes)
  
})

gsea_hallmark <- do.call(rbind, gsea_hallmark)
gsea_hallmark <- gsea_hallmark[!is.na(gsea_hallmark$pval), ]
qsave(gsea_hallmark, file="/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/drug/gsea/sts_hallmark_pathway_drug.qs")
write.csv(gsea_hallmark[, -9], file = "/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/drug/gsea/sts_hallmark_pathway_drug.csv", row.names = FALSE)

###################################################
## GSEA: STS gene association result and KEGG
###################################################

dat <- qread(file="/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/drug/gene_drug_assoc_sts_late_integrated_celligner.qs")
dat <- dat[!is.na(dat$Coef), ]
drug <- unique(dat$drug)

gsea_kegg <- lapply(1:length(drug), function(i){
  
  df <- dat[dat$drug == drug[i], ] 
  #ranks <- df$Coef
  ranks <- sign(df$Coef) * (-log10(df$Pval))
  names(ranks) <-df$Gene
  ranks <- sort(ranks, decreasing = T)
  
  fgseaRes <- fgsea(kegg_pathway, 
                    ranks, 
                    minSize=15, 
                    maxSize = 500)
  fgseaRes <- fgseaRes[!is.na(fgseaRes$padj), ]
  
  data.frame(drug = drug[i], 
             fgseaRes)
  
})

gsea_kegg <- do.call(rbind, gsea_kegg)
gsea_kegg <- gsea_kegg[!is.na(gsea_kegg$pval), ]

qsave(gsea_kegg, file="/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/drug/gsea/sts_kegg_pathway_drug.qs")
write.csv(gsea_kegg[, -9], file = "/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/drug/gsea/sts_kegg_pathway_drug.csv", 
          row.names = FALSE)


