########################################
## load libraries
########################################

library(qs)
library(PharmacoGx)
library(data.table)
library(magrittr)
library(meta)
library(metafor)

##################################################################
## define meta-analysis function
##################################################################
# meta-analysis 
metafun <- function(coef, se, study, pval, n, cancer.type){
  
  data <- data.frame( Gene = feature,
                      Study = as.character( study ),
                      N = n,
                      Coef = as.numeric(as.character( coef )),
                      SE = as.numeric(as.character( se )),
                      Pval = as.numeric(as.character( pval )),
                      Cancer_type = as.character( cancer.type ))
  
  data <- data[ order( data$Coef ) , ]
  
  ## at least 3 studies needed to do the random effect meta-analyses
  if(nrow(data) >= 3){
    
    meta <- metagen( TE = Coef,
                     seTE = SE ,
                     data = data ,
                     studlab = Study ,
                     fixed = FALSE ,
                     random = TRUE ,
                     control = list( maxiter = 10000 , stepadj=0.5 ) )
    
    meta_res <- data.frame(Gene = feature,
                           Coef = meta$TE.random ,
                           SE = meta$seTE.random ,
                           CI_lower = meta$lower.random ,
                           CI_upper = meta$upper.random ,
                           Pval = meta$pval.random ,
                           I2 = meta$I2 ,
                           Q_Pval = meta$pval.Q )
  }else{
    
    print("not enough studies to do meta-analysis")
    meta <- NA
    meta_res <- data.frame(Gene = feature,
                           Coef = NA,
                           SE =  NA,
                           CI_lower = NA,
                           CI_upper = NA,
                           Pval = NA,
                           I2= NA,
                           Q_Pval = NA)
    }
  
  return(list(input_data = data,
              meta_output = meta,
              meta_summery = meta_res))
 
}

##################################################################
## setup directory
##################################################################

dir_data <- '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/data'
dir_drug <- '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/result/drug/'
dir_aligned <- '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/result/aligned/'
dir_output <- '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/result/drug/'

#######################################################
## load gene-drug data (after alignment)
#######################################################

dat <- qread( file.path(dir_drug, "drug_rna_aligned_sts.qs"))
pset_aac <- dat$pset_aac 
pset_aligned <- dat$pset_aligned
gene_ann <- dat$gene_ann

# CCLE cellines
id <- colnames(pset_aac )[grep("-ccle", colnames(pset_aac ))]
pset_aac <- pset_aac[, id]
pset_aligned <- pset_aligned[, id]

# remove drugs with more than 70% missing values across CLs
missing_aac <- sapply(1:nrow(pset_aac), function(k){
  length(which(is.na(pset_aac[k,])))
})

missing_df <- data.frame(perc = round((missing_aac/ncol(pset_aac)) * 100), drug= rownames(pset_aac)) 
missing_df <- missing_df[missing_df$perc < 70, ]
pset_aac <- pset_aac[rownames(pset_aac) %in% missing_df$drug, ] # 5 drugs were removed


gene_drug_res <- lapply(1:nrow(pset_aac), function(k){ 
  
  print(k)
  dat_aac <- pset_aac[k, ]
  
  res <- lapply(1:nrow(pset_aligned), function(i){
    
    fit <- cor.test(as.numeric(pset_aligned[i,]), dat_aac)

    data.frame(Ensembl_ID = gene_ann[gene_ann$Symbol == rownames(pset_aligned)[i], "EnsemblGeneId"],
               gene_name = rownames(pset_aligned)[i], 
               drug = rownames(pset_aac)[k],
               estimate = fit$estimate,
               df = fit$parameter + 2,
               lower = fit$conf.int[1],
               upper = fit$conf.int[2],
               pvalue = fit$p.value)
    
  })  
  
  cor_res <- do.call(rbind, res)
  cor_res$padj <- p.adjust(cor_res$pvalue, method = "BH")
  cor_res
  
})

gene_drug_assoc <- do.call(rbind, gene_drug_res)
rownames(gene_drug_assoc) <- NULL

qsave(gene_drug_assoc, file= file.path(dir_output, "gene_drug_assoc_sts_ccle_ctrp.qs"))

sig_res <- gene_drug_assoc[gene_drug_assoc$padj < 0.15, ]
write.csv(sig_res, file=file.path(dir_output, "gene_drug_assoc_sts_ccle_ctrp_sigFDR.csv"))

#######################################################
## GDSC
#######################################################
pset_aac <- dat$pset_aac 
pset_aligned <- dat$pset_aligned
gene_ann <- dat$gene_ann

# GDSC cellines
id <- colnames(pset_aac )[grep("-gdsc", colnames(pset_aac ))]
pset_aac <- pset_aac[, id]
pset_aligned <- pset_aligned[, id]

# remove drugs with more than 70% missing
missing_aac <- sapply(1:nrow(pset_aac), function(k){
  length(which(is.na(pset_aac[k,])))
})

missing_df <- data.frame(perc = round((missing_aac/ncol(pset_aac)) * 100), drug= rownames(pset_aac)) 
missing_df <- missing_df[missing_df$perc < 70, ]

pset_aac <- pset_aac[rownames(pset_aac) %in% missing_df$drug, ] # 5 drugs were removed

gene_drug_res <- lapply(1:nrow(pset_aac), function(k){ 
  
  print(k)
  dat_aac <- pset_aac[k, ]
  
  res <- lapply(1:nrow(pset_aligned), function(i){
    
    fit <- cor.test(as.numeric(pset_aligned[i,]), dat_aac)
    data.frame(Ensembl_ID = gene_ann[gene_ann$Symbol == rownames(pset_aligned)[i], "EnsemblGeneId"],
               gene_name = rownames(pset_aligned)[i], 
               drug = rownames(pset_aac)[k],
               estimate = fit$estimate,
               df = fit$parameter + 2,
               lower = fit$conf.int[1],
               upper = fit$conf.int[2],
               pvalue = fit$p.value)
    
  })  
  
  cor_res <- do.call(rbind, res)
  cor_res$padj <- p.adjust(cor_res$pvalue, method = "BH")
  cor_res
  
})

gene_drug_assoc <- do.call(rbind, gene_drug_res)
rownames(gene_drug_assoc) <- NULL

qsave(gene_drug_assoc, file=file.path(dir_output, "gene_drug_assoc_sts_gdsc.qs"))

sig_res <- gene_drug_assoc[gene_drug_assoc$padj < 0.15, ]
write.csv(sig_res, file=file.path(dir_output, "gene_drug_assoc_sts_gdsc_sigFDR.csv"))

#######################################################
## NCI-sarcoma
#######################################################

pset_aac <- dat$pset_aac # it needs to be changed
pset_aligned <- dat$pset_aligned
gene_ann <- dat$gene_ann

# NCI-sarcoma cellines
id <- colnames(pset_aac )[grep("-nci", colnames(pset_aac ))]
pset_aac <- pset_aac[, id]
pset_aligned <- pset_aligned[, id]

# remove drugs with more than 70% missing
missing_aac <- sapply(1:nrow(pset_aac), function(k){
  length(which(is.na(pset_aac[k,])))
})

missing_df <- data.frame(perc = round((missing_aac/ncol(pset_aac)) * 100), drug= rownames(pset_aac)) 
missing_df <- missing_df[missing_df$perc < 70, ]

pset_aac <- pset_aac[rownames(pset_aac) %in% missing_df$drug, ] # 5 drugs were removed

gene_drug_res <- lapply(1:nrow(pset_aac), function(k){ 
  
  print(k)
  dat_aac <- pset_aac[k, ]
  
  res <- lapply(1:nrow(pset_aligned), function(i){
    
    fit <- cor.test(as.numeric(pset_aligned[i,]), dat_aac)
    data.frame(Ensembl_ID = gene_ann[gene_ann$Symbol == rownames(pset_aligned)[i], "EnsemblGeneId"],
               gene_name = rownames(pset_aligned)[i], 
               drug = rownames(pset_aac)[k],
               estimate = fit$estimate,
               df = fit$parameter + 2,
               lower = fit$conf.int[1],
               upper = fit$conf.int[2],
               pvalue = fit$p.value)
    
  })  
  
  cor_res <- do.call(rbind, res)
  cor_res$padj <- p.adjust(cor_res$pvalue, method = "BH")
  cor_res
  
})

gene_drug_assoc <- do.call(rbind, gene_drug_res)
rownames(gene_drug_assoc) <- NULL

qsave(gene_drug_assoc, file=file.path(dir_output, "gene_drug_assoc_sts_nci.qs"))

sig_res <- gene_drug_assoc[gene_drug_assoc$padj < 0.15, ]
write.csv(sig_res, file=file.path(dir_output, "gene_drug_assoc_sts_nci_sigFDR.csv"))

####################################################################
## meta-analysis
####################################################################

res.ccle <- qread(file=file.path(dir_output, "gene_drug_assoc_sts_ccle_ctrp.qs"))   
res.gdsc <- qread(file=file.path(dir_output, "gene_drug_assoc_sts_gdsc.qs"))  
res.nci <- qread(file=file.path(dir_output, "gene_drug_assoc_sts_nci.qs"))  

res.ccle$study <- "CCLE/CTRP"
res.gdsc$study <- "GDSC"
res.nci$study <- "NCI"
res <- rbind(res.ccle, res.gdsc, res.nci)

drugs <- intersect(intersect(res.ccle$drug, res.gdsc$drug),
                   res.nci$drug)
genes <- intersect(intersect(res.ccle$gene_name, res.gdsc$gene_name),
                   res.nci$gene_name)

meta.res.drug <- lapply(1:length(drugs), function(k){
  
  print(k)
  df <- res[res$drug == drugs[k], ]
  meta.res <- lapply(1:length(genes), function(i){
    
    df.genes <- df[df$gene_name == genes[i], ]
    
    fit <- metacor(df.genes$estimate, 
                   df.genes$df, 
                   sm = "cor",
                   control=list(stepadj=0.5, maxiter=1000))
    
   data.frame(gene_name = df.genes$gene_name[1],
              Ensembl_ID = df.genes$Ensembl_ID[1],
              r = fit$TE.random,
              se = fit$seTE.random,
              pval = fit$pval.random,
              I2 = fit$I2)
    
  })
  
  meta.res <- do.call(rbind, meta.res)
  meta.res$padj <- p.adjust(meta.res$pval, method = "BH")
  meta.res$drug <- drugs[k]
  meta.res
  
})

meta.res.drug <- do.call(rbind, meta.res.drug)
sig_res <- meta.res.drug[meta.res.drug$padj < 0.15, ]
write.csv(sig_res, file= file.path(dir_output, "gene_drug_assoc_sts_meta_sigFDR.csv"))
qsave(meta.res.drug, file= file.path(dir_output, "gene_drug_assoc_sts_meta.qs"))


