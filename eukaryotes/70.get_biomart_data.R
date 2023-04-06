library("biomaRt")
library('stringr')
source("constant_variables.R")


worker <- function(row, parameters) {
  
  species <- as.character(row$species_abbreviation)
  abbreviation <- as.character(row$abbreviation)
  class <- as.character(row$ensembl_class)
  dataset_name <- as.character(row$dataset)
  
  f <- paste("./tsv/", abbreviation, "_mapping_GO", ".tsv", sep="") 
  #if (file.exists(f)) {
  #  return()
  #}

  host <- as.character(parameters$host)
  mart <- as.character(parameters$mart)
  port <- as.character(parameters$port)

  if (mart == 'ensembl'){
    # new way - use 'useEnsembl' function to define the mirror
    dataset <- useEnsembl(biomart = mart, 
                          dataset = dataset_name, 
                          mirror = "useast")
  } else {
    dataset <- useMart(biomart=mart, dataset=dataset_name, host=host, port=port)
  }


  # TRANSLATION
  if (species %in% specific_species) {
    trans_attrs <- strsplit(str_replace(paste(trans_attrs, collapse= ','), 'ensembl_gene_id', 'external_gene_name'), ',')[[1]]
    #annotation_attrs <- strsplit(str_replace(paste(annotation_attrs, collapse= ','), 'ensembl_gene_id', 'external_gene_name'), ',')[[1]]
    go_attrs <- strsplit(str_replace(paste(go_attrs, collapse= ','), 'ensembl_gene_id', 'external_gene_name'), ',')[[1]]
  }
  
  translation_query <- getBM(attributes=trans_attrs, mart=dataset)
  output_file <- paste("./tsv/", abbreviation, "_translation", ".tsv", sep="")
  write.table(translation_query, file=output_file, sep="\t", na="", quote=FALSE, row.names=FALSE)
  
  # ANNOTATION
  #annotation_query <- getBM(attributes=annotation_attrs, mart=dataset)
  #output_file <- paste("./tsv/", abbreviation, "_annotation", ".tsv", sep="")
  #write.table(annotation_query, file=output_file, sep="\t", na="", quote=FALSE, row.names=FALSE)
  
  # GO
  go_query <- getBM(attributes=go_attrs, mart=dataset)
  output_file <- paste("./tsv/", abbreviation, "_mapping_GO", ".tsv", sep="")
  write.table(go_query, file=output_file, sep="\t", na="", quote=FALSE, row.names=FALSE)
}


###
### DESCRIPTION:
### 
### This script retrieves the necessary data for the subsequent semantic analysis 
### of proteostasis related gene sets, from the Ensembl database. 
###
### For each species, three types of data frames are constructed:
### 
### 1. Translation: Mapping of ensembl ids and gene symbols
### 2. Annotation: Description of genes
### 3. GO: Mapping the GO terms - gene symbols (or ensembl ids)
### 
### Output folder: 'tsv'
###

#annotation_attrs <- c("ensembl_gene_id", "start_position", "end_position", 
#                      "chromosome_name", "band", "strand", "description")
trans_attrs <- c("ensembl_gene_id", "ensembl_gene_id", "ensembl_gene_id")
go_attrs <- c("ensembl_gene_id", "go_id", "go_linkage_type", "namespace_1003")
specific_species <- c('hsapiens', 'mmusculus', 'celegans', 'dmelanogaster', 'drerio', 
                     'ggallus', 'scerevisiae', 'athaliana', 'scerevisiae')

dir.create('./tsv/', showWarnings = FALSE)
df <- read.csv('./files/final_species_set_for_analysis.tsv', sep='\t')

a <- sapply(1:nrow(df), function(i) {
  parameters <- ensembl_parameters_df[ensembl_parameters_df$class == as.character(df[i,]$ensembl_class),]
  worker(df[i,], parameters)
})

