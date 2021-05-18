library("biomaRt")
library('stringr')

worker <- function(row) {

    species <- as.character(row$species)
    class <- as.character(row$class)
    mart <- as.character(row$mart)
    host <- as.character(row$host)
    extension <- as.character(row$dataset_extension)
    dataset <- useMart(mart, dataset=paste(species, extension, sep=""), host=host)
    translation_query <- getBM(attributes=trans_attrs, mart=dataset)
    main_dir <- paste('./data/',class,'/preprocessing/annotation/model_organisms/additional_data/', sep="")
    dir.create(main_dir, showWarnings = FALSE, recursive = TRUE)
    filename <- paste(main_dir,species,'_translation','.tsv', sep="")
    write.table(translation_query, file=filename, sep="\t", na="", quote=FALSE, row.names=FALSE)
    
}


# DESCRIPTION:
# 
# That script calls a request function to send a query to the Ensembl database 
# in order to retrieve alternative gene names for each gene symbol of model
# species. That task is necessary, as the manually curated gene lists of proteostasis
# contain various types of gene names (primary and secondary symbols,
# ensembl ids, entrez ids etc), so all of them should be translated into primary 
# gene symbols.

source("constant_variables.R")
trans_attrs = c("external_gene_name", "ensembl_gene_id")
a <- sapply(as.character(model_species_df$species), function(species){
      row <- model_species_df[model_species_df$species == species,]
      worker(row)
})





