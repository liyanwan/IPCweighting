# Install the TCGAbiolinks package
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("TCGAbiolinks")

# Load the necessary packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(pheatmap)
library(purrr)

# Query for LGG gene expression data
query_lgg <- GDCquery(project = "TCGA-LGG",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      experimental.strategy = 'RNA-Seq',
                      workflow.type = "STAR - Counts",
                      access = 'open')
# Download the data
GDCdownload(query_lgg)
# Prepare the data
lgg_data = GDCprepare(query_lgg)
rownames(lgg_data) = rowData(lgg_data)$gene_name
lgg_data = lgg_data[!duplicated(rownames(lgg_data)),]
data = lgg_data[rowData(lgg_data)$gene_type=="protein_coding", ]
# The columns of gene_expression_matrix represent observations, and the rows represent genes
gene_expression_matrix = assay(data)

sample_type = substr(colnames(gene_expression_matrix), 14, 15)
primary_tumor_samples = colnames(gene_expression_matrix)[sample_type == "01"]
genes_matrix_primary = gene_expression_matrix[, colnames(gene_expression_matrix) %in% primary_tumor_samples]
genes_matrix_t = t(genes_matrix_primary)

# Query for clinical metadata of TCGA-CGG
query_clinical <- GDCquery(project = "TCGA-LGG",
                           data.category = "Clinical",
                           data.type = "Clinical Supplement",
                           data.format = "BCR XML") # Age and sex
GDCdownload(query_clinical)
clinical_data <- GDCprepare_clinic(query_clinical, "patient")
# Case TCGA-R8-A6YH is not included in clinical_data because its file is in data.type = pathology_report
# nrow(clinical_data) == 515 # Should be TRUE
clinical_data$bcr_patient_barcode = as.character(clinical_data$bcr_patient_barcode)
# results = getResults(query_clinical)

sample_barcodes = substr(rownames(genes_matrix_t), 1, 12)
clinical_barcodes = clinical_data$bcr_patient_barcode
common_barcodes = intersect(sample_barcodes, clinical_barcodes) # sort(common_barcodes) == sort(clinical_barcodes) should be TRUE
genes_subset = genes_matrix_t[sample_barcodes%in%common_barcodes, ]
clinical_data_subset = clinical_data[clinical_data$bcr_patient_barcode%in%common_barcodes, ]
genes_subset_barcodes = substr(rownames(genes_subset), 1, 12)
genes_subset_reorder = genes_subset[match(clinical_data_subset$bcr_patient_barcode, genes_subset_barcodes), ]  
all(substr(rownames(genes_subset_reorder), 1, 12) == clinical_data_subset$bcr_patient_barcode) #should be TRUE

time_point = 5
setset = clinical_data_subset
X = data.frame(genes_subset_reorder)
Y_initial = data.frame(bcr_patient_barcode = setset$bcr_patient_barcode,
                       year_of_initial_pathologic_diagnosis = setset$year_of_initial_pathologic_diagnosis, 
                       vital_status = setset$vital_status, 
                       days_to_death = setset$days_to_death, 
                       days_to_last_followup = setset$days_to_last_followup)
Y_initial$vital_status = as.character(Y_initial$vital_status)
# Remove observations that have missing survival times, censor times, and status
# Remove observations that have missing age values
# Remove observations that have negative and 0 observed times
na_info = which(is.na(Y_initial$vital_status) & is.na(Y_initial$days_to_death) & is.na(Y_initial$days_to_last_followup))
if(length(na_info) > 0){
  Y_initial = Y_initial[-na_info, ]
  X = X[-na_info, ]
}
empty_info = which((Y_initial$vital_status == "") & is.na(Y_initial$days_to_death) & is.na(Y_initial$days_to_last_followup))
if(length(empty_info) > 0){
  Y_initial = Y_initial[-empty_info, ]
  X = X[-empty_info, ]
}
# missing_age = which(is.na(X$age))
# if(length(missing_age) > 0){ # Should be FALSE
#   Y_initial = Y_initial[-missing_age, ]
#   X = X[-missing_age, ]
# }
mismatch_dead = which((Y_initial$vital_status == "Dead") & (is.na(Y_initial$days_to_death)) & (is.na(Y_initial$days_to_last_followup)))
if(length(mismatch_dead) > 0){
  Y_initial = Y_initial[-mismatch_dead, ]
  X = X[-mismatch_dead, ]
}
mismatch_alive = which((Y_initial$vital_status == "Alive") & (is.na(Y_initial$days_to_death)) & (is.na(Y_initial$days_to_last_followup)))
if(length(mismatch_alive) > 0){
  Y_initial = Y_initial[-mismatch_alive, ]
  X = X[-mismatch_alive, ]
}
Y = Y_initial
Y$vital_status = ifelse(!is.na(Y$days_to_death), "Dead", Y$vital_status)
Y$observed_time = ifelse(Y$vital_status == "Dead", Y$days_to_death, Y$days_to_last_followup)
Y$observed_time = Y$observed_time/365.25
neg_obs = which(Y$observed_time <= 0) 
if(length(neg_obs) > 0){
  Y = Y[-neg_obs, ]
  X = X[-neg_obs, ]
}
Y$sigma = ifelse(Y$vital_status == "Dead", 1, 0)
Y$E = ifelse(Y$observed_time <= time_point & Y$sigma == 1, 1, ifelse(Y$observed_time > time_point, 0, NA)) # dim(Y) = 467 X 8

# Add a column of IDH with 1p19q Subtype
# lgg.gbm.subtype <- TCGAquery_subtype(tumor = "lgg")
IDH_info = read.csv("IPCweighting/LGG/IDH.csv", stringsAsFactors = FALSE)
X_withPredictors <- as.data.frame(X)
case_names = substr(rownames(X_withPredictors), 1, 12)
X_withPredictors[["IDH/1p19q Subtype"]] = IDH_info$IDH.codel.subtype[match(case_names, IDH_info$Case)] #dim(X_withscaledPredictors) = 467 X 78

all(substr(rownames(X_withPredictors), 1, 12) == Y$bcr_patient_barcode) # Should be true
Y_importantInfo = data.frame(observed_time = Y$observed_time, 
                             sigma = Y$sigma, 
                             E = Y$E)

X_genes <- X_withPredictors[,1:19934]
# Normalize library sizes to account for sequencing depth differences and
# Compute the corresponding log-transformed counts
dge <- DGEList(counts = t(X_genes))
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
log_cpm <- cpm(dge, normalized.lib.sizes = TRUE, log=TRUE)

# Get the top 300 genes with highest variance
gene_var <- apply(log_cpm, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE)[1:300])
top_300_genes_matrix = t(log_cpm[top_genes, ])
top_300_genes_scaled = as.data.frame(scale(top_300_genes_matrix))
genes_with_IDH = data.frame(top_300_genes_scaled, "IDH/1p19q Subtype" = X_withPredictors[["IDH/1p19q Subtype"]])
all(substr(rownames(top_300_genes_scaled), 1, 12) == Y_importantInfo$bcr_patient_barcode)

saveRDS(lgg_data, file = "~/IPCweighting/lgg_data.rds")
saveRDS(clinical_data, file = "~/IPCweighting/clinical_data.rds")
saveRDS(top_300_genes_scaled, file = "~/IPCweighting/genes_top300.rds")
saveRDS(genes_with_IDH, file = "~/IPCweighting/300genes_IDH.rds")
saveRDS(Y, file = "~/IPCweighting/Y_allinfo.rds")
saveRDS(Y_importantInfo, file = "~/IPCweighting/Y.rds")
