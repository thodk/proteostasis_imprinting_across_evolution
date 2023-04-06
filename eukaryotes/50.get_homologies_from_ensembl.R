library("biomaRt")
source('constant_variables.R')


worker <- function(reference_species, parameters) {

    class <- as.character(parameters$class)
    mart <- as.character(parameters$mart)
    suffix <- as.character(parameters$dataset_extension)
    host <- as.character(parameters$host)
    port <- as.character(parameters$port)
    
    main_dir <- define_main_dir(class)
    output_dir <- paste(main_dir,'annotation/homologies/', sep='')
    dir.create(output_dir, showWarnings = FALSE)
    
    species_df <- read.csv(paste(main_dir,'03_species_with_rRNA_and_HSPs_annotation.tsv', sep=''), 
                           sep='\t')

    for (tmp_species in species_df$species_abbreviation) {
        
        print(c(reference_species, tmp_species))
        
        f <- paste(output_dir, tmp_species, '_', reference_species, '_homologs.tsv', sep='')  
        if (file.exists(f)) {
            next
        }
        if (class == 'ensembl' && tmp_species %in% model_species) {
            homology_attrs = c('external_gene_name', 'ensembl_gene_id', 
                               'organism_homolog_associated_gene_name')
        } else if (class == 'ensembl') {
            homology_attrs = c('external_gene_name', 'ensembl_gene_id', 
                               'organism_homolog_ensembl_gene')
        } else if (tmp_species %in% model_species) {
            homology_attrs = c('external_gene_name', 'ensembl_gene_id', 
                               'organism_eg_homolog_associated_gene_name')   
        } else {
            homology_attrs = c('external_gene_name', 'ensembl_gene_id', 
                               'organism_eg_homolog_ensembl_gene')         
        }

        homology_attrs[3] <- gsub('organism', tmp_species, homology_attrs[3])
        dataset_name <- paste(reference_species, suffix, sep="")
        
        if (mart == 'ensembl'){
          # new way - use 'useEnsembl' function to define the mirror
          mart_obj <- useEnsembl(biomart = mart, 
                                dataset = dataset_name, 
                                mirror = "useast")
          #mart_obj <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://asia.ensembl.org",
          #                    dataset=dataset_name)
        } else {
          mart_obj <- useMart(biomart=mart, dataset=dataset_name, host=host, port=port)
        }
        
        
        connection = FALSE
        k = 0
        while (connection == FALSE && k < 10) {
            tryCatch({
                  df <- getBM(mart = mart_obj, attributes = homology_attrs)
                  write.table(x=df, file=f, sep='\t')
                  connection = TRUE
              },  error = function(e) {
                  k <- k + 1
              }, finally = {
                  k <- k + 1
            })
        }
        
    }

}


###
### DESCRIPTION:
### 
### This script retrieves the genomic homologies of model (reference) species
### with all the other species of the same ensembl class (taxonomy). The query 
### is  performed using the appropriate attribute names of each dataset.
###
### The output tsv file contains three columns:
### 1. ref species gene symbol
### 2. ref species ensembl id
### 3. other species gene symbol
### 
### example:
### PLDALPHA4	AT1G55180	AMTR_s00005p00260980
### MGD3	AT2G11810	AMTR_s00033p00201820
###
### Output directory: /data/class/preprocessing/annotation/homologies/
###
###

a <- sapply(as.character(ensembl_parameters_df$class), function(class) {

      rows <- model_species_df[model_species_df$class == class,]
      for (species in as.character(rows$species)) {
        #if (species != 'drerio') {
        #  next
        #}
        row <- model_species_df[model_species_df$species == species,]
        worker(species, row)
      }

})








