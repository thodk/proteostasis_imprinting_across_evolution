library("biomaRt")
library('stringr')


worker <- function(row, parameters) {
  
  species <- as.character(row$species_abbreviation)
  abbreviation <- as.character(row$abbreviation)
  class <- as.character(row$ensembl_class)
  dataset_name <- as.character(row$dataset)
  
  host <- as.character(parameters$host)
  mart <- as.character(parameters$mart)
  port <- as.character(parameters$port)

  dataset <- useMart(biomart=mart, dataset=dataset_name, host=host, port=port)
  
  # TRANSLATION
  if (species %in% specific_species) {
    trans_attrs <- strsplit(str_replace(paste(trans_attrs, collapse= ','), 'ensembl_gene_id', 'external_gene_name'), ',')[[1]]
    annotation_attrs <- strsplit(str_replace(paste(annotation_attrs, collapse= ','), 'ensembl_gene_id', 'external_gene_name'), ',')[[1]]
    go_attrs <- strsplit(str_replace(paste(go_attrs, collapse= ','), 'ensembl_gene_id', 'external_gene_name'), ',')[[1]]
  }
  
  translation_query <- getBM(attributes=trans_attrs, mart=dataset)
  output_file <- paste("./tsv/", abbreviation, "_translation", ".tsv", sep="")
  write.table(translation_query, file=output_file, sep="\t", na="", quote=FALSE, row.names=FALSE)
  
  # ANNOTATION
  annotation_query <- getBM(attributes=annotation_attrs, mart=dataset)
  output_file <- paste("./tsv/", abbreviation, "_annotation", ".tsv", sep="")
  write.table(annotation_query, file=output_file, sep="\t", na="", quote=FALSE, row.names=FALSE)
  
  # GO
  go_query <- getBM(attributes=go_attrs, mart=dataset)
  output_file <- paste("./tsv/", abbreviation, "_mapping_GO", ".tsv", sep="")
  write.table(go_query, file=output_file, sep="\t", na="", quote=FALSE, row.names=FALSE)
}



# DESCRIPTION:
# 
# That script retrieves the necessary data for the subsequent analysis from the 
# Ensembl database. Three types of data frames:
# 
# 1. Translation: Mapping of ensembl ids and gene symbols
# 2. Annotation: Description of genes
# 3. GO: contains the GO terms - gene symbols (or ensembl ids) mapping
# Output directory: /tsv/

source("constant_variables.R")
annotation_attrs <- c("ensembl_gene_id", "start_position", "end_position", 
                      "chromosome_name", "band", "strand", "description")
trans_attrs <- c("ensembl_gene_id", "ensembl_gene_id", "ensembl_gene_id")
go_attrs <- c("ensembl_gene_id", "go_id", "go_linkage_type", "namespace_1003")
specific_species <- c('hsapiens', 'mmusculus', 'celegans', 'dmelanogaster', 'drerio', 
                     'ggallus', 'scerevisiae', 'athaliana', 'scerevisiae', 'spombe')

dir.create('./tsv/', showWarnings = FALSE)
df <- read.csv('./files/final_species_set_for_analysis.tsv', sep='\t')

a <- sapply(1:nrow(df), function(i) {
  parameters <- parameters_df[parameters_df$class == as.character(df[i,]$ensembl_class_class),]
  worker(df[i,], parameters)
  1
})

