##########################################################################################
###############   FUNCTIONS CREATED TO RETRIEVE, EDIT AND EXPORT CaSQ      ###############
###############   Simple Interaction Files (SIF) FROM CellDesigner TO      ###############
###############   HiPathia "ATT" and "SIF" FORMAT USING Minerva AS PROXY   ###############
##########################################################################################


if (!require(pacman)) {
  install.packages("pacman")
}

if (!require(BiocManager)) {
  install.packages("BiocManager")
}

if (!requireNamespace("AnnotationDbi", quietly = TRUE)){
  BiocManager::install("AnnotationDbi")
}
  
if (!requireNamespace("hipathia", quietly = TRUE)){
  BiocManager::install("hipathia", version = "3.10")
}

if (!requireNamespace("SummarizedExperiment", quietly = TRUE)){
  BiocManager::install("SummarizedExperiment", version = "3.10")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)){
  BiocManager::install("org.Hs.eg.db", version = "3.10")
}


pacman::p_load("here","igraph","rentrez","reutils","hipathia", "biomaRt", "utils", "stringr", "SummarizedExperiment", "magrittr",
               "AnnotationDbi","org.Hs.eg.db","dplyr","tidyr", "openxlsx", "data.table", "gdata", "httr", "jsonlite")


##################################################
## Project: COVID-19 Disease Map
## Script purpose: Translate raw CellDesigner SIF obtained from CaSQ to Entrez identifiers using MINERVA 
## Date: 11.06.2020
## Author: Marek Ostaszewski
##################################################

## A convenience function to handle API queries
ask_GET <- function(furl, fask) {
  print(URLencode(paste0(furl, fask)))
  resp <- httr::GET(url = URLencode(paste0(furl, fask)), write_memory(),
                    httr::add_headers('Content-Type' = "application/x-www-form-urlencoded"),
                    ### Currently ignoring SSL!
                    httr::set_config(config(ssl_verifypeer = 0L)))
  if(httr::status_code(resp) == 200) {
    ### when the content is sent as a zip file, it needs to be handled differently,
    ### i.e. saved to a tmp file and unzipped
    if(headers(resp)$`content-type` == "application/zip") {
      tmp <- tempfile()
      tmpc <- file(tmp, "wb")
      writeBin(httr::content(resp, as = "raw"), con = tmpc)
      close(tmpc)
      unzipped <- unzip(tmp)
      file.remove(tmp)
      return(unzipped)
    }
    return(httr::content(resp, as = "text"))
  }
  return(NULL)
}

### An utility function to retrieve Entrez based on the species id
### if the id is a complex, the function goes recursively and fetches the ids of elements in this complex
group_elements <- function(feid, felements, fentrez) {
  pos <- which(felements$elementId == feid)
  ### Any elements that may be nested in the 'feid' (CellDesigner alias)
  incs <- felements$elementId[felements$complexId %in% felements$id[pos]]
  if(length(incs) > 0) {
    ### If nested elements found, run the function recursively for the contained elements
    return(paste(unlist(sapply(incs, group_elements, felements, fentrez)), collapse = ";"))
  } else {
    ### If no nested elements, return Entrez
    rid <- fentrez[[feid]]
    if(is.na(rid)) {
      ### If Entrez not available, return name
      rid <- felements$name[pos]
    }
    return(rid)
  }
}

### A workaround function to get information about hypothetical complexes; 
### currently MINERVA API does not support this, we need to get the entire CD file and parse it
get_groups <- function(fname) {
  message(paste0("Getting groups for ", fname, "..."))
  library(xml2)
  ## Currently comment out the MINERVA download, some content dlded as binary!
  cd_map <- read_xml(ask_GET(mnv_base,
                             paste0("models/",
                                    models$idObject[models$name == fname],
                                    ":downloadModel?handlerClass=lcsb.mapviewer.converter.model.celldesigner.CellDesignerXmlParser")))
  ### CellDesigner namespace
  ns_cd <- xml_ns_rename(xml_ns(read_xml("<root>
                                            <sbml xmlns = \"http://www.sbml.org/sbml/level2/version4\"/>
                                            <cd xmlns = \"http://www.sbml.org/2001/ns/celldesigner\"/>
                                          </root>")), 
                         d1 = "sbml", d2 = "cd")
  ### Get complex ids
  cids <- xml_attr(xml_find_all(cd_map, "//cd:complexSpeciesAlias", ns_cd), "species")
  ### For each check, which is hypothetical
  hypocs <- sapply(cids, 
                   function(x) {
                     tcid <- xml_find_first(cd_map, paste0("//sbml:species[@id='", x, "']"), ns_cd)
                     ifelse(length(tcid) > 0, xml_text(xml_find_first(tcid, ".//cd:hypothetical", ns_cd)), NA)
                   })
  names(hypocs) <- gsub("s_id_", "",names(hypocs))
  return(hypocs)
}

### A convenience function that makes use of the previous ones to retrieve the translated SIF file from CellDesigner
### and CaSQ using Minerva. The function allows to apply the operation to all the diagrams listed in "sif_raw_gitLab" 
#' @param sif_raw_gitLab list of raw sif files obtained from CaSQ.

extract_translatedSIF_Minerva <- function(sif_raw_gitLab, 
                                          base_url = "https://git-r3lab.uni.lu/covid/models/-/raw/master/Executable%20Modules/SBML_qual_build/sif/"){
  
  message(paste0("Translating..", sif_raw_gitLab))
  
  ### Read in the raw SIF version (here straight from the gitLab Covid19)
  
  raw_sif <- read.table(url(paste0(base_url, sif_raw_gitLab)),
                        sep = " ", header = F, stringsAsFactors = F)
  
  diagram <- gsub(pattern = "_raw.sif", ".xml", sif_raw_gitLab)
  
  
  diag_name <- res[str_detect(res$Resource, diagram), "Name"]
  
  ### Get elements of the chosen diagram
  model_elements <- fromJSON(ask_GET(paste0(mnv_base,"models/",models$idObject[models$name == diag_name],"/"), 
                                     "bioEntities/elements/?columns=id,name,type,references,elementId,complexId,bounds"),
                             flatten = F)
  
  message("Fetching entrez ids...")
  
  
  ### Get information about Entrez identifiers from MINERVA elements
  entrez <- sapply(model_elements$references, function(x) ifelse(length(x) == 0, NA, x[x$type == "ENTREZ", "resource"]))
  names(entrez) <- model_elements$elementId
  
  cgroups <- get_groups(diag_name)
  
  message("Translating...")
  
  
  ### Create a copy
  translated_sif <- raw_sif
  
  
  ### Retrieve Entrez and type for the entire columns of sources and targets
  s.entrez <- sapply(raw_sif[,1], group_elements, model_elements, entrez)
  s.type <- sapply(raw_sif[,1], function(x) { ifelse(x %in% names(cgroups),
                                                     ifelse(is.na(cgroups[x]), "complex", "group"),
                                                     "node") })
  t.entrez <- sapply(raw_sif[,3], group_elements, model_elements, entrez)
  t.type <- sapply(raw_sif[,3], function(x) { ifelse(x %in% names(cgroups),
                                                     ifelse(is.na(cgroups[x]), "complex", "group"),
                                                     "node") })
  ### Collect x.y information
  s.xy <- t(sapply(raw_sif[,1], function(x) unlist(model_elements$bounds[model_elements$elementId == x, c(3,4)])))
  colnames(s.xy) <- c("source.x", "source.y")
  t.xy <- t(sapply(raw_sif[,3], function(x) unlist(model_elements[model_elements$elementId == x,1][,c(3,4)])))
  colnames(t.xy) <- c("targets.x", "targets.y")
  
  ### Combine into a single data frame
  translated_sif <- data.frame(source = s.entrez, source.type = s.type,
                               sign = raw_sif[,2], 
                               target = t.entrez, target.type = t.type,
                               s.xy, t.xy, stringsAsFactors = F)
  
  ## In case export is desired
  # write.table(translated_sif, file = paste0(diag_name,"_translated_sif.txt"),
  #             sep = "\t", quote = F, col.names = T, row.names = F)
  
  return(translated_sif)
  message("Done.")  
  
}


#################################################
## Project: COVID-19 Disease Map
## Script purpose: Translate Entrez identifiers SIF  files from MINERVA to HIPATHIA SIF and ATT files for metagraphinfo compilation.
## Date: 21.10.2020 
## Author: Marina Esteban Medina
##################################################


###################################################################################################
### Function to pass translatedSIF_Minerva sif files information to Hipathia att and sif files ####
#' @param minerva_translated_sif the sif file from CaSQ with entrez identifiers from MINERVA
#' @param diagram_name name of the diagram we are adding
#' @param outputfolder the output file path, it will return the sif an att file in Hipathia format
#' @param exceptions file of the nodes exceptions, such as viral genes with entrez ID, that Hipathia does not processes
########################################################################################################################

hipathia.att.sif <- function (minerva_translated_sif, diagram_name, outputfolder, 
                              exceptions = file.path(paste0(base_dir, "/Resources/Hipathia/sif_att_converter/resources/exceptions.txt"))) {
  
  message(paste0("Working on...", diagram_name))

  exceptions_2avoid <- read.table(exceptions, stringsAsFactors = F) %>% mutate_all(as.character) %>% .$V1
    
  ## Clean auto interaction of nodes
  minerva_translated_sif <- mutate_all(minerva_translated_sif, as.character)
  minerva_translated_sif <- minerva_translated_sif[!minerva_translated_sif$source == minerva_translated_sif$target, ]
  
  ## Clean nodes that are not genes with and entrez ID nor a Function node.
  ## We clean the source and target columns to avoid non gene-gene interactions.
  minerva_translated_sif  <- minerva_translated_sif[which(!grepl("[a-z];[0-9]+|[A-Z];[0-9]+|[0-9]+;[a-z]|[0-9]+;[A-Z]|[A-Z]+[a-z]+[0-9]+[[:punct:]];[0-9]+|Unfolded protein", minerva_translated_sif$source)),] %>% 
    .[which(!grepl("[a-z];[0-9]+|[A-Z];[0-9]+|[0-9]+;[a-z]|[0-9]+;[A-Z]|[A-Z]+[a-z]+[0-9]+[[:punct:]];[0-9]+|Unfolded protein", .$target)),] 
  
  ## Replace space with "_" in the diagrams name
  
  diagram_name <- gsub(" ", "_", diagram_name)
  
  ## Start bulding the att data frame, cols = ID, LABEL, X, Y, COLOR, SHAPE, TYPE, LABEL.CEX, LABEL.COLOR, WIDTH, HEIGHT, genesList
  nodes <- unique(c(minerva_translated_sif$source, minerva_translated_sif$target))
  
  ## Create a table with the type of nodes, complex or protein family
  
  nodes_gene_type <- data.frame(nodes = c(minerva_translated_sif$source, minerva_translated_sif$target),
                                type = c(minerva_translated_sif$source.type, minerva_translated_sif$target.type))%>% 
    .[!duplicated(.$nodes),] %>% mutate_all(as.character)
  
  
  
  id_num <- seq_along(nodes)
  
  ## Create empty vectors for the columns of the att data.frame that we will then export to ".att" file of each diagram
  ID = sapply(id_num, function(x) {paste("N", paste0("hsaCovid19", diagram_name ), x , sep = "-")})
  
  sif <- data.frame(ID, nodes)
  
  label <- c()
  shape <- c()
  type <- c()
  genesList <- c()
  
  for (i in seq_along(nodes)){
    
    ## In the case of human gene nodes with entrez identifiers
    if (grepl ("[A-Z]|[a-z]", nodes[i]) == F && nodes[i] != "8673700" && !(nodes[i] %in% exceptions_2avoid) ){ ## those are new entrez ID of SARS-CoV2 which we want to avoid
      
      message(paste0("Translating...", nodes[i]))
      
      label[i] <- str_split(nodes[i], pattern = ";")  %>% sapply ( . , function(x) { 
        mapIds(org.Hs.eg.db, keys = x, column = "SYMBOL", keytype = "ENTREZID")
      }) %>% paste(., collapse = ";")
      
      shape[i] <- "ellipse"
      
      type[i] <- "gene"
      
      genesList[i] <- nodes[i]
      
      if(nodes_gene_type$type[nodes_gene_type$nodes %in% nodes[i]] == "complex"){
        genesList[i] <- gsub( ";", ",/,", nodes [i]) ##  substitute group sign (;) by what Hipathia takes as complex sign (,/,)
      }
      
      ## In the case of function nodes 
      
    }else if(grepl ("^[A-Z]+[aeiou]+[[:punct:]]*|^[a-z]+[aeiou]+[[:punct:]]*", nodes[i]) == T){ 
      
      label[i] <- nodes[i]
      
      shape[i] <- "rectangle" ## Hipathia symbol for "function" element in the graph
      
      type[i] <- "node"
      
      genesList[i] <- "NA"
      
      ID[i] <- paste0(ID[i],"_func")
      
      
      
    }else{
      
      message(paste0("Translating...", nodes[i]))
      
      ## We only use this if we want to keep the non gene nodes visually  
      label[i] <- nodes[i]
      
      shape[i] <- "circle"
      
      type[i] <- "other"
      
      genesList[i] <- "NA"
    }
  }
  
  
  
  X <- c()
  
  Y <- c()
  
  for (i in seq_along(nodes)){
    
    if(nodes[i] %in% minerva_translated_sif$source){
      
      X[i] <- minerva_translated_sif$source.x[match(nodes[i], minerva_translated_sif$source)] 
      
      Y[i] <- minerva_translated_sif$source.y[match(nodes[i], minerva_translated_sif$source)] 
      
    }else{
      
      X[i] <- minerva_translated_sif$targets.x[match(nodes[i], minerva_translated_sif$target)] 
      
      Y[i] <- minerva_translated_sif$targets.y[match(nodes[i], minerva_translated_sif$target)] 
      
    }
    
  }
  
  if(length(id_num) == 0 ){
    
    message(paste0("No human genes founded in ", diagram_name, "...DELETED "))
    
  }else{
    
    ## Color, label.cex, label.color, width and height are constant for all the nodes
    att_hipathia = data.frame ( ID, label, X, Y, color = "white", shape, type, label.cex = "0.62", label.color = "black", 
                                width = "46", height = "17", genesList )
    
    
    ## Edit the .sif file in Hipathia format, subsitute interactions for "activation" or "inhibition".
    
    interaction = gsub("POSITIVE", "activation", minerva_translated_sif$sign) %>% gsub("NEGATIVE", "inhibition", .)
    
    sif_hipathia <- data.frame(source = sif$ID [match(minerva_translated_sif$source, sif$nodes)], interaction,
                               target = sif$ID [match(minerva_translated_sif$target, sif$nodes)])
    
    sif_hipathia <- sif_hipathia[sif_hipathia$source %in% att_hipathia$ID, ]
    
    ## Create the "name.pathways_hsa.txt" file. This file serves as list of the pathways included in the metaginfo (metagraph-database) Hipathia takes as pathways database.
    
    name.pathways_hsa_file <- data.frame(code = paste0("Covid19", diagram_name), name = diagram_name)
    
    
    ## Export the "att", "sif" and  "name.pathways_hsa.txt"  to the output folder selected
    
    write.table(att_hipathia, file = paste(outputfolder, paste0(paste0("hsaCovid19", diagram_name ) ,".att"), sep = "/"), sep = "\t", quote = F, col.names = T, row.names = F )
    
    write.table(sif_hipathia, file = paste(outputfolder, paste0(paste0("hsaCovid19", diagram_name ) ,".sif"), sep = "/"), sep = "\t", quote = F, col.names = F, row.names = F )
    
    write.table(name.pathways_hsa_file, file = paste(outputfolder, "name.pathways_hsa.txt", sep = "/"),append = T, sep = "\t", col.names = F, row.names = F, quote = F )
    
    message(paste0(diagram_name, "... Done"))
    
    ## The function also returns the  "sif and "att" data.frames created in a list for downstream actions.
    
    return(list("att" = att_hipathia, "sif" = sif_hipathia))
    
    
  }
  
  
}

