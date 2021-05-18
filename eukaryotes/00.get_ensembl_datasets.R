library("biomaRt")
library("stringr")


worker <- function(row) {
    output_dir = paste('./data/',row$class,'/preprocessing/files/', sep='')
    dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
    print(as.character(row$mart))

    mart <- useMart(biomart=as.character(row$mart), host=as.character(row$host), 
                    port=as.integer(row$port))

    datasets_df <- listDatasets(mart)
    datasets_df$abbreviation <- apply(datasets_df, 1, function(row) {
        unlist(strsplit(row['dataset'], split='_'))[1]
    })
    
    write.table(x=datasets_df, 
		file=paste(output_dir, 'ensembl_datasets.tsv', sep=''),
		sep='\t')
}



source('constant_variables.R')

sapply(as.character(ensembl_parameters_df$class), function(c) {
        row <- ensembl_parameters_df[ensembl_parameters_df$class == c,]
        worker(row)
})




