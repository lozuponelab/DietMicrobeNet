# load necessary packages 
packages <- c("readr", "dplyr", "argparse", "tidyr")
for (package in packages) {
  if(!require(package, character.only = T)){
    install.packages(package)
    library(package)
  }
}

# create parser
parser <- ArgumentParser()
parser$add_argument('--diet_file', help = "CSV file containing food items")
parser$add_argument('--content_file', help = 'CSV file containing content info on food')
parser$add_argument('--ExDes_file', help = 'CSV file containing external descriptors on compounds')
parser$add_argument('--meta_o_file', help = 'file path for compound metadata output')
args <- parser$parse_args()

#' From foods found in FooDB get list of KEGG compounds found in their collective metabolomes
#'
#' @param diet_df CSV generated through included website 
#' @param content_df content CSV  from FooDB
#' @param external_descriptor_df CompoundExternalDescriptor CSV from FooDB
#' @param meta_out_path output path for metadata file on compounds
#'
#' @return file with a list of KEGG Compounds 
#' 
get_diet_fooDB_compounds <- function(diet_df, 
                                     content_df, 
                                     external_descriptor_df, 
                                     meta_out_path){ 
  
  # load in data 
  foods <- read_csv(diet_df)
  Content <- read_csv(content_df)
  CExtDes <- read_csv(external_descriptor_df)
  
  # get food IDs
  ids <- foods$id

  # Get compound IDs, I believe these are labeled as the Source IDs in the contents file
  food_content <- Content[Content$food_id %in% ids, ]
  food_content <- food_content %>% 
    select(food_id, orig_food_common_name, source_id)
  
  comp_ids <- unique(food_content$source_id)

  # get external descriptors
  extDes <- CExtDes[CExtDes$compound_id %in% comp_ids, ]
  extDes <- extDes %>%
    select(external_id, compound_id) %>%
    na.omit()

  # get only external descriptors that are KEGG related
  comp_KEGG <- unique(subset(extDes, grepl("^C[0-9]{5}$", external_id)))
  
  # create a metadata about the compounds 
  colnames(foods) <- c("food_id", "name", "name_scientific", "food_frequency")
  meta <- food_content %>% 
    select(source_id, food_id, orig_food_common_name) %>% 
    rename_with(~c('compound_id', 'food_id', 'common_name'),
                c(source_id, food_id, orig_food_common_name)) %>% 
    left_join(comp_KEGG, by = 'compound_id') %>% 
    left_join(foods, by = 'food_id') %>% 
    drop_na(external_id)
  
  colnames(meta) <- c("foodb_compound_id", "food_id", "org_food", "kegg_id", "name", "name_scientific", "food_frequency")
  
  write.csv(meta, file = meta_out_path, row.names = F)
}

# call function
get_diet_fooDB_compounds(diet_df = args$diet_file,
                        content_df = args$content_file,
                        external_descriptor_df = args$ExDes_file,
                        meta_out_path = args$meta_o_file)

