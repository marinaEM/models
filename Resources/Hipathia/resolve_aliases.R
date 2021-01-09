##################################################
## Project: COVID-19 Disease Map
## Script purpose: Translate raw CellDesigner SIF to Entrez identifiers using MINERVA 
## Date: 05.06.2020
## Author: Marek Ostaszewski
##################################################

library(httr)
library(jsonlite)

### A convenience function to handle API queries
ask_GET <- function(furl, fask, unpack = T) {
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

### Define the source file (GitLab, raw link)
diagram <- "https://git-r3lab.uni.lu/covid/models/-/raw/master/Curation/ER%20Stress/ER_Stress_stable.xml"

### Read in the raw SIF version (here straight from the github of Aurelien)
raw_sif <- read.table(url("https://git-r3lab.uni.lu/covid/models/-/raw/master/Executable%20Modules/SBML_qual_build/sif/ER_Stress_stable_raw.sif"),
                      sep = " ", header = F, stringsAsFactors = F)

### Read the list of resources to be integrated, from the MINERVA build scripts
res <- read.csv(url("https://git-r3lab.uni.lu/covid/models/raw/master/Integration/MINERVA_build/resources.csv"),
                header = T, stringsAsFactors = F)

diag_name <- res[res$Resource == diagram, "Name"]

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
models <- ask_GET(mnv_base, "models/")
models <- fromJSON(models, flatten = F)

this_refs <- models[models$name == diag_name]

### Get elements of the chosen diagram
model_elements <- fromJSON(ask_GET(paste0(mnv_base,"models/",models$idObject[models$name == diag_name],"/"), 
                                   "bioEntities/elements/?columns=id,name,type,references,elementId,complexId,bounds"),
                           flatten = F)

message("Fetching entrez ids...")
### Get information about Entrez identifiers from MINERVA elements
entrez <- sapply(model_elements$references, function(x) ifelse(length(x) == 0, NA, x[x$type == "ENTREZ", "resource"]))
names(entrez) <- model_elements$elementId

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
                             s.xy, t.xy)

write.table(translated_sif, file = "translated_sif.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)
message("Done.")
