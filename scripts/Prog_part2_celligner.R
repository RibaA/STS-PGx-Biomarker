################################################################################
## libraries
################################################################################

source("/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/scripts/Prog_celligner_function.R")

##################################################################
## setup directory
##################################################################

dir_input <- '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/data'
dir_output <- '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/result/aligned/'
dir_annot <- '/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/aligned/'

#################################################################
## run Celligner
#################################################################
# Download gene annotation (ensemble ids) 

hgnc.complete.set <- data.table::fread( file.path(dir_annot, "hgnc_complete_set.txt")) %>% as.data.frame() 

# load corrected cell line and tumor patient data
pgx_rna <- run_Celligner(data_dir =  dir_annot, 
                         dat =  file.path(dir_input, "PGx_corrected_gse_rna_sts.qs"), 
                         remove_cPCA_dims = 1:7, 
                         plotname = " ")

qsave(pgx_rna, file = file.path(dir_output, "pgx_rna_celligner_sts.qs"))

################################################################################
## visualization Part 1: before alignment
################################################################################

plot_uncorrected_data_type(org_dat = qread(file.path(dir_input, "PGx_corrected_gse_rna_sts.qs")), 
                           before_plot = "pgx_tumor_uncorrected",
                           dir_output = dir_output)

plot_uncorrected_data_study(org_dat = qread(file.path(dir_input, "PGx_corrected_gse_rna_sts.qs")), 
                            before_plot = "pgx_tumor_uncorrected_study",
                            dir_output = dir_output)

################################################################################
## visualization Part 2: after alignment
################################################################################

dat <- qread(file.path(dir_output, "pgx_rna_celligner_sts.qs"))

Celligner_alignment_plot_type(aligned_data = dat$comb_obj, 
                              after_plot = "pgx_tumor_corrected",
                              dir_output = dir_output)

Celligner_alignment_plot_study(aligned_data = dat$comb_obj, 
                               org_dat = qread(file.path(dir_input, "PGx_corrected_gse_rna_sts.qs")),
                               after_plot =  "pgx_tumor_corrected_study")



