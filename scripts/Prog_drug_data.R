########################################
## load libraries
########################################

library(qs)
library(PharmacoGx)
library(VennDiagram)

##################################################################
## setup directory
##################################################################

dir_data <- '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/data'
dir_aligned <- '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/result/aligned/'
dir_output <- '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/result/drug/'

################################################################################
## load drug response data
################################################################################
dat <- qread(file= file.path(dir_data, "PGx_sarc_data.qs"))

gene_ann <- dat$rna$CCLE$marray_gene_ann

ctrp_aac <- dat$aac$CTRP$drug_aac
gdsc_aac <- dat$aac$GDSCv2$drug_aac
nci_aac <- dat$aac$NCI$drug_aac


## Venn-diagram

n1 <- length(rownames(ctrp_aac))
n2 <- length(rownames(gdsc_aac))
n3 <- length(rownames(nci_aac))
n12 <- length(intersect(rownames(ctrp_aac), rownames(gdsc_aac)))
n13 <- length(intersect(rownames(ctrp_aac), rownames(nci_aac)))
n23 <- length(intersect(rownames(gdsc_aac), rownames(nci_aac)))
n123 <- length(intersect(intersect(rownames(ctrp_aac), rownames(gdsc_aac)), rownames(nci_aac)))

pdf(file = file.path(dir_output, "venn_drug.pdf"), width = 6, height = 6)

venn.plot <- draw.triple.venn(
  area1 = n1,
  area2 = n2,
  area3 = n3,
  n12 = n12,
  n13 = n13,
  n23 = n23,
  n123 = n123,
  category = c("CTRP", "GDSC", "NCISarcoma"),
  fill = c("#64894DFF", "#99B6BDFF", "#ECC9A0FF"),
  cat.col = c("#252525", "#252525", "#252525"),
  #lty = "dashed",
  cex = 2,
  cat.cex = 1.5
)

dev.off()


# filter to common drugs across three PSets
int <- intersect(intersect(rownames(ctrp_aac), rownames(gdsc_aac)),
                 rownames(nci_aac)) # 80 drugs

ctrp_aac <- ctrp_aac[int, ]
colnames(ctrp_aac) <- paste(colnames(ctrp_aac), "ccle", sep="-")
nci_aac <- nci_aac[int, ]
colnames(nci_aac) <- paste(colnames(nci_aac), "nci", sep="-")
gdsc_aac <- gdsc_aac[int, ]
colnames(gdsc_aac) <- paste(colnames(gdsc_aac), "gdsc", sep="-")

pset_aac <- cbind(ctrp_aac, gdsc_aac, nci_aac)

## load aligned data 
dat <- qread(file.path(dir_aligned, "pgx_rna_celligner_sts.qs"))

dat_aligned <- as.data.frame(Seurat::GetAssayData(dat$comb_obj))
dat_ann <- data.frame(sampleID = dat$comb_obj$sampleID,
                      type = dat$comb_obj$type)

cl_id <- dat_ann[dat_ann$type == "CL", "sampleID"]
pset_aligned <- dat_aligned[ , cl_id]

int <- intersect(colnames(pset_aligned), colnames(pset_aac)) # 61 cell lines in common
pset_aac <- pset_aac[, colnames(pset_aac) %in% int] # 28 cell lines (include three drugs) and 24 (gdsc and ctrp)
pset_aligned <- pset_aligned[, colnames(pset_aligned) %in% int]

dat <- list(pset_aac = pset_aac, pset_aligned = pset_aligned, gene_ann = gene_ann)
qsave(dat, file= file.path(dir_output, "drug_rna_aligned_sts.qs"))

