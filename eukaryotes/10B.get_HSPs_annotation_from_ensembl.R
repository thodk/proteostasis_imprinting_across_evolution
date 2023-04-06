library("biomaRt")
library("stringr")
source('constant_variables.R')

str_modification <- function(s) {
    list_of_letters <- lapply(unlist(strsplit(s, '')), function(x) {
        if (suppressWarnings(!is.na(as.numeric(x))) || x==' '){
            x
        } else {
            paste(c('[',x,toupper(x),']'), collapse='') 
        }
    })
    paste(unlist(list_of_letters), collapse='')
}


create_df <- function(df, keywords) {
    keywords <- unique(unlist(lapply(keywords , function(s) {str_modification(s)})))
    s <- paste0(keywords, collapse="|")
    prohibited_s <- paste0(prohibited_keywords, collapse="|")
    hsp_df <- df[which(str_detect(as.vector(df$description), s ) == TRUE), ]
    final_df <- hsp_df[which(str_detect(as.vector(hsp_df$description), prohibited_s) == FALSE), ]
    return(final_df)
}


worker <- function(row) {
  
    main_dir <- define_main_dir(as.character(row$class))
    hsp40_dir = paste(main_dir,'ensembl_data/hsp40_raw_files/', sep='')
    hsp70_dir = paste(main_dir,'ensembl_data/hsp70_raw_files/', sep='')
    dir.create(paste(main_dir,'ensembl_data/', sep=''), showWarnings = FALSE)    
    dir.create(hsp40_dir, showWarnings = FALSE)
    dir.create(hsp70_dir, showWarnings = FALSE)
    
    datasets_df <- read.csv(paste(main_dir,'files/','ensembl_datasets.tsv', sep=''), sep='\t')
    datasets <- datasets_df$dataset

    for (dataset in datasets) {
        species <- unlist(strsplit(dataset, '_'))[1]

      	# try sequential attempts to connect to the server of ensembl (max 10) 
      	# because sometimes the connections cannot be established 
        mart = as.character(row$mart)
        host = as.character(row$host)
        port = as.integer(row$port)
        connection = FALSE
        k = 0
        while (connection == FALSE && k < 10){
            print(c(species, k))
            tryCatch({
                dataset <- useMart(biomart=mart, dataset=dataset, host=host, port=port)
                connection = TRUE
                k <- k + 1               
            }, error = function(e) {
                k <- k + 1
            },
            finally={
                k <- k + 1
            }
          )
        }
        
      	# 1. the derived data frame should contain an ensembl_id and the respective
      	# gene symbol ('external_gene_name') for each entry.
      	# 2. the chunk below examines if the list of dataset attributes contains 
      	# the field 'external_gene_name' else it uses the 'ensembl_gene_id'
        tryCatch({
            attrs <- c('ensembl_gene_id', 'external_gene_name',  'description')
            df <- getBM(attrs, mart=dataset)
        }, error = function(e) {
            attrs <- c('ensembl_gene_id', 'ensembl_gene_id',  'description')
            df <- getBM(attrs, mart=dataset)
            names(df) <- c('ensembl_gene_id', 'external_gene_name',  'description')
        })

      	# 1. check if the request has returned data and store it in tsv files, 
        # labeled as 'raw'
        # 2. load keywords and filter out entries that do not fulfil the requirements
        if (class(df) == "function"){
            next
        } else {
            final_df <- create_df(df, valid_keywords_for_hsp70)
            print(dim(final_df))            
            f <- paste(hsp70_dir, species, '.tsv', sep='')
            write.table(final_df, file=f, sep='\t', quote=FALSE)
            final_df <- create_df(df, valid_keywords_for_hsp40)
            f <- paste(hsp40_dir, species, '.tsv', sep='')
            write.table(final_df, file=f, sep='\t', quote=FALSE)
	      }
    }
}


###
### DESCRIPTION:
### 
### This script retrieves the annotation for HSPs from the Ensembl database. 
### The worker function is executed for the species of each ensembl class 
### separately. The final goal is to detect which species have adequate 
### annotation for HSPs and filter out all the others. The output lists of 
### annotated HSPs will be combined with similar data retrieved from the 
### UniProt database (see 10A.get_HSPs_annotation_from_uniprot.py).
###
###

a <- sapply(as.character(ensembl_parameters_df$class), function(c){
  	row <- ensembl_parameters_df[ensembl_parameters_df$class == c,]
  	worker(row)
})





