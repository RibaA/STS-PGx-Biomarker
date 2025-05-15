########################################
## load libraries
########################################

library(qs)
library(PharmacoGx)
library(data.table)
library(magrittr)
library(corrplot)
library(ggplot2)

##################################################################
## define function
##################################################################

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

#########################################################################################################################
#########################################################################################################################
########################################## Load Gene Drug response (AAC) Association ####################################
#########################################################################################################################
#########################################################################################################################
## class of drugs
load("/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Data/GDSC.RData")
dat_drug <- drugInfo(dat)

## univariable results: STS 
dat <- qread("/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/drug/gene_drug_assoc_sts_late_integrated_celligner_lm_7pcs.qs")
dat_drug <- dat_drug[rownames(dat_drug) %in% unique(dat$drug), ]

## get the class of drugs
dat_drug$TARGET_PATHWAY <- ifelse(dat_drug$TARGET_PATHWAY %in% c("Other, kinases", "RTK signaling///Other, kinases", "Other"), 
                                  "Other", dat_drug$TARGET_PATHWAY)

dat_drug$TARGET_PATHWAY <- ifelse(dat_drug$TARGET_PATHWAY %in% c("PI3K/MTOR signaling"), 
                                  "PI3K-MTOR signaling", dat_drug$TARGET_PATHWAY)
dat_drug$TARGET_PATHWAY <- ifelse(dat_drug$TARGET_PATHWAY %in% c("RTK signaling///IGF1R signaling"), 
                                  "RTK signaling-IGF1R signaling", dat_drug$TARGET_PATHWAY)
dat_drug$TARGET_PATHWAY <- ifelse(dat_drug$TARGET_PATHWAY %in% c("RTK signaling///EGFR signalin"), 
                                  "RTK signaling-EGFR signalin", dat_drug$TARGET_PATHWAY)


target_pathway <- unique(dat_drug$TARGET_PATHWAY)
freq_target_pathway <- lapply(1:length(target_pathway), function(k){
  
 data.frame(target_pathway = target_pathway[k], 
            freq = nrow( dat_drug[dat_drug$TARGET_PATHWAY == target_pathway[k], ] ))
  
})

freq_target_pathway  <- do.call(rbind, freq_target_pathway)
freq_target_pathway <- freq_target_pathway[order(freq_target_pathway$freq, decreasing = TRUE), ]
freq_target_pathway_included <- freq_target_pathway[freq_target_pathway$freq >= 2, ]
target_pathway_included <- freq_target_pathway_included$target_pathway

dat_drug_class <- lapply(1:length(target_pathway_included), function(k){
  
  dat_drug[dat_drug$TARGET_PATHWAY == target_pathway_included[k], "treatmentid"]
  
})

names(dat_drug_class) <- target_pathway_included

#########################################################################################################################
#########################################################################################################################
#########################################  Pearson similarity metric  ###################################################
#########################################################################################################################
#########################################################################################################################

cor_res <- lapply(1:length(dat_drug_class), function(k){
  
  sub_res <- dat[dat$drug %in% dat_drug_class[[k]], ]
  
  sub_res_matrix <- lapply(1:length(dat_drug_class[[k]]), function(j){
    
    sub_res[sub_res$drug == dat_drug_class[[k]][j], "Coef"]
    
  })
  
  sub_res_matrix <- do.call(cbind, sub_res_matrix)
  colnames(sub_res_matrix) <- dat_drug_class[[k]]
  res <- rcorr(as.matrix(sub_res_matrix))
  data.frame(target_pathway = names(dat_drug_class)[k],
             flattenCorrMatrix(res$r, res$P) )
  
  
})

cor_res <- do.call(rbind, cor_res)

write.csv(cor_res, file="/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/drug/cor_pcl_lm.csv")

## box plot for all the pharmacological classes

cor_res <- read.csv("/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/drug/cor_pcl_lm.csv")
#cor_res$target_pathway <- ifelse(cor_res$target_pathway %in% c("Other, kinases", "RTK signaling///Other, kinases", "Other"), "Other", cor_res$target_pathway)

#cor_res$label <- ifelse(cor_res$target_pathway %in% c("Other, kinases", "Other"), 0, 1)
#cor_res <- cor_res[cor_res$label ==1, ]

jpeg(file="/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/drug/boxplot_lm.jpg", 
     width = 500, height = 500, res = 150)

p <- ggplot(cor_res, aes(x=reorder(target_pathway, -cor), y=cor)) + 
  geom_boxplot(width = 0.5, fill= "#F4E7C5FF", color = "#78847FFF") +
  coord_flip() + 
  ylim(c(-0.5, 1)) +
  xlab("") +
  ylab("Pearson correlation \n (predictive value)") +
  theme(axis.text.x=element_text(size=8),
        axis.title=element_text(size=9),
        axis.text.y=element_text(size=8),
        strip.text = element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none",
        legend.text = element_text(size = 6, face="bold"),
        legend.title = element_blank()) 

p

dev.off()

## upper triangles

for(k in 1:length(dat_drug_class)){
  
  print(k)
  
  sub_res <- dat[dat$drug %in% dat_drug_class[[k]], ]
  
  sub_res_matrix <- lapply(1:length(dat_drug_class[[k]]), function(j){
    
    sub_res[sub_res$drug == dat_drug_class[[k]][j], "Coef"]
    
  })
  
  sub_res_matrix <- do.call(cbind, sub_res_matrix)
  rem <- which(dat_drug_class[[k]] == "Unii-40E3azg1MX")
  if(length(rem) > 0){
    
    dat_drug_class[[k]][rem] <- "BMS-536924"
    
  }
  colnames(sub_res_matrix) <- dat_drug_class[[k]]
  # for RTK pathway
  
  #sub_res_matrix <- sub_res_matrix[, colnames(sub_res_matrix) %in% c("Axitinib", "Linifanib", "Quizartinib")]
  
  res <- cor(sub_res_matrix)
  
 # if(k == 2){ names(dat_drug_class)[k] <- "PI3K_MTOR signaling"   } 
 #  if(k == 5){ names(dat_drug_class)[k] <- "Other_kinases"   }  
 # if(k == 13){ names(dat_drug_class)[k] <- "RTK signaling Other kinases"   } 
  
  jpeg(file=paste("/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/drug/PCL/",
                  paste(names(dat_drug_class)[k], ".jpg", sep=""), sep=""),
       width = 650, height = 650, res=150)
  
  p <- corrplot(res, type = "upper", order = "hclust", 
           tl.col = "black", tl.srt = 45, tl.cex = 0.8)
  p
  
  dev.off()
  
}

##########################################################
## all 80 drugs
##########################################################

dat <- qread("/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/drug/gene_drug_assoc_sts_late_integrated_celligner_lm_7pcs.qs")
sig <- dat[dat$padj < 0.05, ]
dat_drug <- dat_drug[rownames(dat_drug) %in% unique(sig$drug), ]

drugs <- unique(sig$drug)
res <- lapply(1:length(drugs), function(k){
  
  df <- dat[dat$drug == drugs[k], "Coef"]
  #df <- ifelse(df > 1, 1, df)
  #df <- ifelse(df < (-1), -1, df)
  df
  
})

res <- do.call(cbind, res)
colnames(res) <- drugs

cor_res <- cor(res)

jpeg(file=paste("/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/drug/PCL/",
                paste("cor_all_lm", ".jpg", sep=""), sep=""),
     width = 1300, height = 1300, res=150)

p <- corrplot(cor_res, type = "upper", order = "hclust", 
              tl.col = "black", tl.srt = 90, tl.cex = 0.8)
p

dev.off()


