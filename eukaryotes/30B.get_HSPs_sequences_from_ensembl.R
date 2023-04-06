library("biomaRt")
library("stringr")
source('constant_variables.R')


collect_data <- function(main_dir, row, species, hsp) {

    # load the file with HSPs annotation from filtered_hsp_annotation
    # and keep the entries which contain the 'ensembl' in the 'resource' column
    df <- read.csv(paste(main_dir,'filtered_HSPs_annotation/',species,'.tsv', sep=''), sep='\t')
    hsp_df <- df[df$hsp == hsp,]
    hsp_df <- hsp_df[hsp_df$resource == 'ensembl',]
    ensembl_gene_ids <- unique(hsp_df$gene_symbol)
    gene_symbols <- unique(hsp_df$gene_symbol)

    # connect to the ensembl server (max attempts = 10)
    mart = as.character(row$mart)
    host = as.character(row$host)
    port = as.integer(row$port)
    
    dataset <- useMart(biomart=mart,
                       dataset=paste(c(species,as.character(row$dataset_extension)), collapse=''),
                       host=host, port=port)

    # the output data frame 
    output_df <- data.frame(gene_symbol=character(), peptide=character())
    
    # firstly try to get sequences for ensembl_gene_id
    if (length(ensembl_gene_ids) > 0) {
          attrs <- c('ensembl_gene_id',  'peptide')
          filters = c('ensembl_gene_id')
          values = ensembl_gene_ids
          connection = FALSE
          k = 0
          while (connection == FALSE && k < 10){
                tryCatch({
                      df1 <- getBM(attrs, values=values, filters=filters, mart=dataset)
                      names(df1) <- c('gene_symbol', 'peptide')
                      output_df <- rbind(output_df, df1)
                      connection = TRUE
                }, error = function(e) {
                      k <- k + 1
                }, finally = {
                      k <- k + 1
                })
          }
    }
  
    # secondly try to get sequences for gene symbols
    if (length(gene_symbols) > 0) {
          attrs <- c('external_gene_name',  'peptide')
          filters <- c('external_gene_name')
          values <- gene_symbols
          connection = FALSE
          k = 0
          while (connection == FALSE && k < 10){
              tryCatch({
                    df2 <- getBM(attrs, values=values, filters=filters, mart=dataset)
                    names(df2) <- c('gene_symbol', 'peptide')
                    output_df <- rbind(output_df, df2)
                    connection = TRUE
              }, error = function(e) {
                    k <- k + 1
              }, finally = {
                    k <- k + 1
              })
          }
    }

    return(output_df)
}




worker <- function(row) {
    
    main_dir <- define_main_dir(row$class)
    hsp40_dir <- paste(main_dir,'ensembl_data/hsp40_sequences/', sep='')
    hsp70_dir <- paste(main_dir, 'ensembl_data/hsp70_sequences/', sep='')
    dir.create(paste(main_dir,'ensembl_data/', sep=''), showWarnings = FALSE)   
    dir.create(hsp40_dir, showWarnings = FALSE)
    dir.create(hsp70_dir, showWarnings = FALSE)

    # get the species names that have passed the previous steps
    # in order to retrieve their HSPs sequences from ensembl
    species_df <- read.csv(paste(main_dir,'03_species_with_rRNA_and_HSPs_annotation.tsv', sep=''), 
                             sep='\t', row.names='X')

    for (r in 1:nrow(species_df)) {
        species <- as.character(species_df[r, 'species_abbreviation'])
        abbreviation <- as.character(species_df[r, 'abbreviation'])
        f1 <- paste(hsp40_dir, abbreviation, '.tsv', sep='')
        f2 <- paste(hsp70_dir, abbreviation, '.tsv', sep='')
        
        # retrieve HSP40 sequences, using the 'collect_data' function
        #if (file.exists(f1) & file.exists(f2)){
        #    next
        #}
        #if (file.exists(f1)) {
            #
        #} else {
            df = collect_data(main_dir, row, species, 'hsp40')
            if (class(df) == "data.frame") {
            		# save the data as tsv and not fasta files
            		# use as file name the specific abbreviation code of each species,
            		# e.g. plant_0, plant_10, fungi_8 etc
                write.table(df, file=f1, sep='\t', quote=FALSE)
            } 
        #}
    
        # retrieve HSP70 sequences, using the 'collect_data' function
        #if (file.exists(f2)) {
            #
        #} else {        
            df = collect_data(main_dir, row, species, 'hsp70')
            if (class(df) == "data.frame") {
            		# save the data as tsv and not fasta files
            		# use as file name the specific abbreviation code of each species,
            		# e.g. plant_0, plant_10, fungi_8 etc
                write.table(df, file=f2, sep='\t', quote=FALSE)
            } 
        #}
    }
}


###
### DESCRIPTION:
### 
### This script retrieves the HSP sequences from the Ensembl database. It reads
### the species-specific files in 'filtered_hsp_annotation' subfolder to get the 
### ids of HSPs sequences and sends a query for them in the Ensembl database. The 
### output files are in tsv format, containing the protein id, gene symbol and 
### sequence string for each entry.
### 
### All the output files are stored in:
### '/data/class/preprocessing//ensembl_data/hsp40_sequences/' &
### '/data/class/preprocessing//ensembl_data/hsp70_sequences/'
###
### Here the species are named with a unique id, based on their ensembl_class, 
### in order to avoid to use their abbreviations, during the analysis. This 
### nomenclature is simpler and provides a more safe way to create the Mongo 
### collections in the following steps:
### 
### example: ensembl_0, ensembl_10, ensembl_20, plant_2, plant_13 etc 
###
###
### Help: this script depends on the connection with the Ensembl repository.
### Potential problems with that connection will kill the process. Usually, the 
### connection is lost after a couple of queries. So, the script should be called
### many times, until its execution without raised connection errors.
###
### Error in useMart(biomart = mart, mirror = "asia", dataset = paste(c(species,  : 
### unused argument (mirror = "asia")
### Calls: sapply -> lapply -> FUN -> worker -> collect_data
### Execution halted
###


a <- sapply(as.character(ensembl_parameters_df$class), function(c){
    row <- ensembl_parameters_df[ensembl_parameters_df$class == c,]
    worker(row)
})






