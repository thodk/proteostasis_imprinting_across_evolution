###
### DESCRIPTION:
###
### This script contains some variables and functions which are used in many 
### different steps of the data acquisition process. Thus, they are concentrated 
### here to loaded thought the command source('constant_variables.R') wherever
### they are necessary. Some of the following variables are related with
### tsv files stored in the 'files' directory.
###
###

ensembl_parameters_df <- read.csv('./files/ensembl_parameters.tsv', sep='\t')
model_species_df <- read.csv('./files/model_species.tsv', sep='\t')
model_species <- as.character(model_species_df$species)

define_main_dir <- function(x) {
    paste('./data/',x,'/preprocessing/', sep='')
}

valid_keywords_for_hsp40 <- c('heat shock 40', 'heat shock protein 40', 
                              'heat shock cognate 40', 'shock 40kDa', 
                              'shock 40 kDa', 'chaperone 40', 'hsp40', 'dnaj')
valid_keywords_for_hsp70 <- c('heat shock 70', 'heat shock protein 70', 
                              'heat shock cognate 70', 'heat shock cognate 71',
                              'shock 70kDa', 'shock 70 kDa', 'shock 71 kDa',
                              'chaperone 70', 'hsp70', 'dnak')
prohibited_keywords <- c('interact', 'cooporate', 'stimulate', 'resemble', 
                         'contain', 'unknown function', 'with the hsp70', 
                         'with the hsp40', 'escort', 'fragment', 'pseudogene', 
                         'complex', 'sense overlapping', 'sense intronic', 
                         'antisense', 'binding', 'exchange factor',
                         'co-chaperone', 'suppressor', 'cooperates with',
                         'organizing protein', 'readthrough', 'divergent transcript')