      ###############################################################################
##########   EXAMPLE CODE: From CellDseigner CaSQ "Raw sif" files  to         ############
##############             Hipathia "sif" and "att" format using Minerva   #############
     #############################################################################
     
###############################################################################################################     
## IMPORTANT: CHANGE THE "base_dir" parameter to the local directory where the repository was cloned !!!!  ####
      
      base_dir <- "~/Desktop/models-master/"
      
#####################################################################################################      
      
#### 1. LOAD LIBRARIES  AND FUNCTIONS  #####

             ## Load functions
      source(file.path(paste0(base_dir, "/Resources/Hipathia/sif_att_converter/functions_retrieve_sif_att_hipathia_from_CellDesigner.R")), echo=TRUE) 

      
#### 2. LOAD FOLDER, CREATE LIST OF RAW_SIF_FILES; LOAD MINERVA MAPS AND INFO NEEDED ####
      
      ### Read the raw sif files from the corresponding folder in the repository
  
      sif_files_path <- file.path(paste0(base_dir,"/Executable Modules/SBML_qual_build/sif"))
      
      raw_sif_list <- list.files(path = sif_files_path , all.files = F) %>% .[grep("*_raw*", .)] 
      
      diagrams <- gsub("_stable_raw.sif", "", raw_sif_list) 
      
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
      
      ### Read the list of resources to be integrated, from the MINERVA build scripts
      res <- read.csv(url("https://git-r3lab.uni.lu/covid/models/raw/master/Integration/MINERVA_build/resources.csv"),
                      header = T, stringsAsFactors = F)
      
      
#### 3. RETRIEVE TTRANSALTED SIF FILES WITH DETAILED INFORMATION OF THE NODES AND INTERACTIONS OF EACH GRAPH USING MINERVA #####
      
  ### The obtained translated_sif_files counts with: X-Y coords info, genesList in entrez Identifier and type of 
  ###  node ( simple, complex or protein family). The info will be retrieved and stored as a data.frame for each diagram.                                                    
      
      translated_sif_files <- lapply(raw_sif_list[-c(6)],
                                     function(x){extract_translatedSIF_Minerva(x)}) 
      
      
#### 4. FORMAT TRANSALTED SIF FILES OBTAINED FOR HIPATHIA METAGRAPHINFO DATABASE ##### 
      
      
      att_sif_Hformat <- mapply(hipathia.att.sif, minerva_translated_sif = translated_sif_files, diagram_name = diagrams[-6], 
                                outputfolder = "~/Desktop/models-master/Resources/Hipathia/sif_att_converter/test")
      
      