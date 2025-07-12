#02_Normalize_Log_Transform_Counts.R
#Normalizes and log transforms RNA seq counts from GTEx v8 tissues

library(DESeq2)
source("path_utils.R")

#define a function to normalize and log transform the raw RNAseq counts
Normalize_log_transform_data <- function(data, drops = c("id", "Name", "Description"), tissue_name = "Tissue") {
  #set the rownames to the ENSG id 
  rownames(data) <- data$Name
  
  #remove all columns except gene counts, safe to do now that the gene names are the rownames
  raw_counts <- data[, !(names(data) %in% drops)]
  
  #Create metadata file
  metadata <- data.frame(colnames(raw_counts))
  metadata$Sample <- tissue_name
  colnames(metadata) <- c("Donor", "Tissue")
  rownames(metadata) <- metadata$Donor
  
  #Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  
  #Get normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  
  #Add a pseudocount of one
  normalized_counts_plus_pseudocount <- normalized_counts + 1
  
  #Log transform the data
  normalized_counts_plus_pseudocount_log_transformed <- log2(normalized_counts_plus_pseudocount)
  
  #Return the log-transformed data
  return(normalized_counts_plus_pseudocount_log_transformed)
}

#Set GTEx raw counts directory
Raw_Counts_path <- file.path(data_dir, "GTEx_Raw_Tissue_RNA_counts")

#list each file to read
Raw_count_files <- list.files(path = Raw_Counts_path, pattern = "\\.gct$", full.names = TRUE)

#Read all count files as a dataframe into a list
raw_counts_list <- lapply(Raw_count_files, function(f) {
  read.csv(f, sep = "\t", header = TRUE)
})

#name each dataframe in the list according to its tissue's name
tissue_names <- sub("gene_reads_2017-06-05_v8_(.*)\\.gct$", "\\1", basename(Raw_count_files))
names(raw_counts_list) <- tissue_names


#normalize and log transform the counts for every tissue
normalized_counts_list <- lapply(raw_counts_list, Normalize_log_transform_data)

#Set GTEx normalized counts directory to export normalized counts
Normalized_Counts_path <- file.path(data_dir, "GTEx_Normalized_Tissue_RNA_counts")

#Export the normalized and log-transformed counts
for (tissue_name in names(normalized_counts_list)) {
  df <- normalized_counts_list[[tissue_name]]
  output_file <- file.path(Normalized_Counts_path, paste0(tissue_name, ".csv"))
  write.csv(df, file = output_file, row.names = TRUE)
}

