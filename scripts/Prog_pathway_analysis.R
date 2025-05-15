########################################
## load libraries
########################################

library(qs)
library(data.table)
library(magrittr)
library(fgsea)
library(reshape2)
#library(enrichR)

##################################################################
## setup directory
##################################################################

dir <- '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/result/drug/'
dir_pathway <- '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/data/'
dir_output <- '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/result/pathway/'

###################################################################
## load KEGG, GO and HALLMARK pathways downloaded from msigdb
##################################################################

hallmark_pathway <- gmtPathways(file.path(dir_pathway, "h.all.v2024.1.Hs.symbols.gmt"))
kegg_pathway <- gmtPathways(file.path(dir_pathway, "c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt"))
go_pathway <- gmtPathways(file.path(dir_pathway, "c5.go.bp.v2024.1.Hs.symbols.gmt"))

dat <- list(KEGG = kegg_pathway, HALLMARK = hallmark_pathway, GO = go_pathway)
qsave(dat, file=  file.path(dir_pathway, "pathway_data.qs"))

###################################################
## GSEA: STS gene association result and HALLMARK
###################################################

dat <- qread(file= file.path(dir, "gene_drug_assoc_sts_meta.qs"))
dat <- dat[!is.na(dat$r), ]
drug <- unique(dat$drug)

gsea_hallmark <- lapply(1:length(drug), function(i){
  
  df <- dat[dat$drug == drug[i], ] 
  ranks <- df$r
  #ranks <- sign(df$r) * (-log10(df$pval))
  names(ranks) <- df$gene_name
  ranks <- sort(ranks, decreasing = T)
  
  fgseaRes <- fgsea(hallmark_pathway, 
                    ranks, 
                    minSize=15, 
                    maxSize = 500)
  
  fgseaRes <- fgseaRes[!is.na(fgseaRes$padj), ]
  fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, paste, collapse = ",")
  data.frame(drug = drug[i], fgseaRes)
  
})

gsea_hallmark <- do.call(rbind, gsea_hallmark)
gsea_hallmark <- gsea_hallmark[!is.na(gsea_hallmark$pval), ]

qsave(gsea_hallmark, file= file.path(dir_output, "hallmark_gsea_pathway_drug.qs"))
write.csv(gsea_hallmark, file = file.path(dir_output, "hallmark_gsea_pathway_drug.csv"), row.names = FALSE)

###################################################
## GSEA: STS gene association result and KEGG
###################################################

gsea_kegg <- lapply(1:length(drug), function(i){
  
  df <- dat[dat$drug == drug[i], ] 
  ranks <- df$r
  #ranks <- sign(df$r) * (-log10(df$pval))
  names(ranks) <- df$gene_name
  ranks <- sort(ranks, decreasing = T)
  
  fgseaRes <- fgsea(kegg_pathway, 
                    ranks, 
                    minSize=15, 
                    maxSize = 500)
  
  fgseaRes <- fgseaRes[!is.na(fgseaRes$padj), ]
  fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, paste, collapse = ",")
  data.frame(drug = drug[i], fgseaRes)
  
})

gsea_kegg <- do.call(rbind, gsea_kegg)
gsea_kegg <- gsea_kegg[!is.na(gsea_kegg$pval), ]

qsave(gsea_kegg, file= file.path(dir_output, "kegg_gsea_pathway_drug.qs"))
write.csv(gsea_kegg, file = file.path(dir_output, "kegg_gsea_pathway_drug.csv"), row.names = FALSE)

###################################################
## GSEA: STS gene association result and GO
###################################################

gsea_go <- lapply(1:length(drug), function(i){
  
  df <- dat[dat$drug == drug[i], ] 
  ranks <- df$r
  #ranks <- sign(df$r) * (-log10(df$pval))
  names(ranks) <- df$gene_name
  ranks <- sort(ranks, decreasing = T)
  
  fgseaRes <- fgsea(go_pathway, 
                    ranks, 
                    minSize=15, 
                    maxSize = 500)
  
  fgseaRes <- fgseaRes[!is.na(fgseaRes$padj), ]
  fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, paste, collapse = ",")
  data.frame(drug = drug[i], fgseaRes)
  
})

gsea_go <- do.call(rbind, gsea_go)
gsea_go <- gsea_go[!is.na(gsea_go$pval), ]

qsave(gsea_go, file= file.path(dir_output, "go_gsea_pathway_drug.qs"))
write.csv(gsea_go, file = file.path(dir_output, "go_gsea_pathway_drug.csv"), row.names = FALSE)

################################################################################
## ORA: STS gene association result and HALLMARK
################################################################################
dat <- qread(file= file.path(dir, "gene_drug_assoc_sts_meta.qs"))
dat <- dat[!is.na(dat$r), ]
drug <- unique(dat$drug)

ora_hallmark <- lapply(1:length(drug), function(i){
  
  df <- dat[dat$drug == drug[i], ] 
  universe <- df$gene_name 
  genes <- df[df$padj < 0.1, "gene_name"]
  
  foraRes <- fora(hallmark_pathway, 
                   genes = genes,
                   universe = universe)
  
  foraRes <- foraRes[!is.na(foraRes$padj), ]
  foraRes$overlapGenes <- sapply(foraRes$overlapGenes, paste, collapse = ",")
  data.frame(drug = drug[i], foraRes)
  
})

ora_hallmark <- do.call(rbind, ora_hallmark)
ora_hallmark <- ora_hallmark[!is.na(ora_hallmark$pval), ]

qsave(ora_hallmark, file= file.path(dir_output, "ora_hallmark_pathway_drug.qs"))
write.csv(ora_hallmark, file = file.path(dir_output, "ora_hallmark_pathway_drug.csv"), row.names = FALSE)


################################################################################
## ORA: STS gene association result and KEGG
################################################################################

ora_kegg <- lapply(1:length(drug), function(i){
  
  df <- dat[dat$drug == drug[i], ] 
  universe <- df$gene_name 
  genes <- df[df$padj < 0.1, "gene_name"]
  
  foraRes <- fora(kegg_pathway, 
                  genes = genes,
                  universe = universe)
  
  foraRes <- foraRes[!is.na(foraRes$padj), ]
  foraRes$overlapGenes <- sapply(foraRes$overlapGenes, paste, collapse = ",")
  data.frame(drug = drug[i], foraRes)
  
})

ora_kegg <- do.call(rbind, ora_kegg)
ora_kegg <- ora_kegg[!is.na(ora_kegg$pval), ]

qsave(ora_kegg, file= file.path(dir_output, "ora_kegg_pathway_drug.qs"))
write.csv(ora_kegg, file = file.path(dir_output, "ora_kegg_pathway_drug.csv"), row.names = FALSE)

################################################################################
## ORA: STS gene association result and GO
################################################################################

ora_go <- lapply(1:length(drug), function(i){
  
  df <- dat[dat$drug == drug[i], ] 
  universe <- df$gene_name 
  genes <- df[df$padj < 0.1, "gene_name"]
  
  foraRes <- fora(go_pathway, 
                  genes = genes,
                  universe = universe)
  
  foraRes <- foraRes[!is.na(foraRes$padj), ]
  foraRes$overlapGenes <- sapply(foraRes$overlapGenes, paste, collapse = ",")
  data.frame(drug = drug[i], foraRes)
  
})

ora_go <- do.call(rbind, ora_go)
ora_go <- ora_go[!is.na(ora_go$pval), ]

qsave(ora_go, file= file.path(dir_output, "ora_go_pathway_drug.qs"))
write.csv(ora_go, file = file.path(dir_output, "ora_go_pathway_drug.csv"), row.names = FALSE)


