# Converter of COVID-19 CellDesigner diagrams "raw SIF" files obtained using CaSQ and MINERVA to HiPathia "SIF" and "ATT" file format.

The present directory contains the code and files necessary to convert *Simple Interaction Files (SIF)*, generated with CaSQ python package, from curated COVID-19 related CellDesigner diagrams to HiPathia's *SIF* and *ATT (attribute)* files. The *ATT* files are generated using Minerva, as a proxy between CaSQ and CellDesigner, to retrieve details of each node. These *SIF* and *ATT* files are needed by the HiPathia mechanistic modelling metagraphinfo database in order to run the algorithm over the diagrams. After running the conversor, each diagram will count with a *"diagram.sif"* and a *"diagram.att"* file,moreover a *"names.pathways.txt"* file listing the names of the diagrams converted will be also generated.

**The directory is structured into 2 folders and 2 scripts.**

### Scripts

** - functions_retrieve_sif_att_hipathia_from_CellDesigner.R**:

The present script contains the functions needed for converting *raw SIF *files from CellDesigner to HiPathia *SIF* and *ATT* files. The script is autodescriptive and counts with 5 functions for the retrieval, preprocessing and formatting of the raw SIF files. The functions enable to handle API querys, retrieve detailed information of each node from a CellDesigner diagram and format the information into HiPathia's SIF and ATT files. Each function is fully describe in the script. DO NOT EDIT (only edit for coding purposes).




** - example_sif_Hipathia.R**

The present script contains the example code to run over the existing diagrams on the [https://git-r3lab.uni.lu/covid/models/-/tree/master/](https://git-r3lab.uni.lu/covid/models/-/tree/master/) repository, please in order to run the code correctly clone the repository on a local directory and recode the "base_dir" variable to that directory. The script will read the *Simple Interaction Files (SIF)* from the CellDesigner diagrams located in */models-master/Executable Modules/SBML_qual_build/sif* in the repository, generate HiPathia formatted *SIF* and *ATT* files for each diagram and a *"names.pathways.txt"* file listing the names of the diagrams converted. The files are by default stored in the **test folder**, this parameter can be change in the *outputfolder* parameter of the *hipathia.att.sif* function.

<u>**CAUTION**:</u> Please when running the script make sure you recode the variable "base_dir" to the local directory where the repo *models-master* was cloned.



### Resources folder
Contains the *"exceptions.txt"* file, which contains for the nodes we desire to avoid, nodes that are not allowed in the metagraphinfo database such as viral gene nodes.

### Test folder
Empty folder that will store the results of running the *"example_sif_Hipathia.R"*


### Contributors

[Marek Ostaszewski](https://fairdomhub.org/people/665)

[Marina Esteban-Medina](https://fairdomhub.org/people/1873)



