##################################################
## Project: COVID-19 Disease Map
## Script purpose: Integrate OmniPathDB with COVID-19 Disease Map on MINERVA 
## Date: 28.05.2020
## Author: Marek Ostaszewski
##################################################

### A convenience function to handle API queries
ask_GET <- function(furl, fask) {
  resp <- httr::GET(url = paste0(furl, fask),
                    httr::add_headers('Content-Type' = "application/x-www-form-urlencoded"))
  if(httr::status_code(resp) == 200) {
    return(httr::content(resp, as = "text"))
  }
  return(NULL)
}

### Load OmniPath
if(!'OmnipathR' %in% installed.packages()[,"Package"]){
  require(devtools)
  install_github('saezlab/OmnipathR')
}

require(OmnipathR)
require(dplyr)

message(paste0("Asking for OmniPath interactions"))

### Get OmniPath interactions
ia_omnipath <- import_Omnipath_Interactions() %>% as_tibble()

### Get MINERVA elements
### The address of the COVID-19 Disease Map in MINERVA
map <- "https://covid19map.elixir-luxembourg.org/minerva/api/"
### Get configuration of the COVID-19 Disease Map, to obtain the latest (default) version
cfg <- fromJSON(ask_GET(map, "configuration/"))
project_id <- cfg$options[cfg$options$type == "DEFAULT_MAP","value"]
### The address of the latest (default) build 
mnv_base <- paste0(map,"projects/",project_id,"/")

message(paste0("Asking for diagrams in: ", mnv_base, "models/"))

### Get diagrams
models <- fromJSON(ask_GET(mnv_base, "models/"), flatten = F)
### Get elements of diagrams
model_elements <- lapply(models$idObject, 
                         function(x) fromJSON(ask_GET(paste0(mnv_base,"models/",x,"/"), 
                                                      "bioEntities/elements/?columns=id,name,type,references,complexId"), 
                                              flatten = F))
names(model_elements) <- models$name

### A convenience function toextract a given annotation type from references
extract_annos <- function(refs, type) {
  refs <- refs[sapply(refs, length) > 0]
  if(length(refs) == 0) { return(NULL) }
  ret <- sapply(refs, function(x) x[x$type == type, "resource"])
  ret <- ret[sapply(ret, function(x) ifelse(is.character(x) & length(x) > 0, TRUE, FALSE))]
}


### A convenience function to parse the annotations
extract_hgnc <- function(model) {
  ### For erroneous response, return NULL
  if(is.null(model)) { return(NULL) }
  ### Only elements that have annotations
  these_refs <- model$references[sapply(model$references, length) > 0]
  ### For empty list, return NULL
  if(length(these_refs) == 0) { return(NULL) }
  ### Get HGNC symbols, for elements that have annotations
  hgncs <- sapply(these_refs, function(x) x[x$type == "HGNC_SYMBOL", "resource"])
  hgncs <- hgncs[sapply(hgncs, function(x) ifelse(is.character(x) & length(x) > 0, TRUE, FALSE))]
  return(unique(unlist(hgncs)))
}

hgncs_of_models <- sapply(model_elements, extract_hgnc)

### Parse HGNCs list and compare to the interactions from the OmniPathDB
for(hmi in 1:length(hgncs_of_models)) {
  hmis <- hgncs_of_models[[hmi]]
  
  message(paste0("At least one interactor in ", names(hgncs_of_models)[hmi]))
  ia_omnipath %>% filter(source_genesymbol %in% hmis | target_genesymbol %in% hmis) -> op_oneint
  print(op_oneint)
  
  message(paste0("Both interactors in ", names(hgncs_of_models)[hmi]))
  ia_omnipath %>% filter(source_genesymbol %in% hmis & target_genesymbol %in% hmis) -> op_twoint
  print(op_twoint)
}

columns <- "id,reactionId,type,products,reactants,modifiers"
model_reactions <- lapply(models$idObject, 
                          function(x) fromJSON(ask_GET(paste0(mnv_base,"models/",x,"/"), 
                                                       paste0("bioEntities/reactions/?columns=", columns)),
                                               flatten = F))


for(mr in 1:length(model_reactions)) {
  print(paste0("Model: ", models$name[mr]))
  ### Get reactant ids (internal MINERVA)
  reas <- sapply(model_reactions[[mr]]$reactants, "[[", "aliasId")
  ### Get modifier ids (internal MINERVA)
  mods <- sapply(model_reactions[[mr]]$modifiers, "[[", "aliasId")
  ### Get product ids (internal MINERVA)
  pros <- sapply(model_reactions[[mr]]$products, "[[", "aliasId")
  if(length(reas) == 0) { next }
  res <- c()
  ### For each reaction
  for(i in 1:length(reas)) {
    ### Get participants
    parts <- unique(c(reas[[i]], mods[[i]], pros[[i]]))
    ### Get HGNC symbols
    parts_hgncs <- extract_annos(model_elements[[mr]][model_elements[[mr]]$id %in% parts, "references"], "HGNC_SYMBOL")
    ### Get OmniPath interactions for these HGNC symbols
    ia_omnipath %>% filter(source_genesymbol %in% parts_hgncs & target_genesymbol %in% parts_hgncs) -> r_match
    if(nrow(r_match) > 0) {
      ### Aggregate results
      res <- rbind(res, r_match)
    }
  }
  print(res)
}
