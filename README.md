## Table of contents
* [General info](#general-info)
* [Usage](#setup)
* [Data Information](#Data-Information)

## General Info: ReactionLumping
Reaction Lumping procedure for metabolic networks for application with thermodynamic metabolic flux analysis (TMFA).
Provided procedure is implemented in MatLab and applicable for genome-scale models.  The procedure is based on identification of groups of metabolites with unknwon \Delta G_f,
such that their eliminations from the model, via rection lumping, can be determined on every group independently. The combined procedure exists out of 2 componentes.
The so called group component tries to eliminate the entire group. If is unsuccesfull, the second component (called sequential lumping) tries every metabolite
individually. For more information please see publication.

## SetUp 
### Use the lumping Function
The function is called 'lumpReactions.m' and can be find in the program folder. Please note that the function is dependent on 'basicLump'. The function
requires a genome-scale model (downloadable for example from the BIGG-database) and an indices list of metabolites with unknown information for metabolites
### Reproduce results 
To reproduce the results described within the publication, as well as an working example please use provided 'main.m' function in which you can
choose between the models provided (via uncommenting).

## Data Information
The redHuman model was retrieved in 10/2020 under: 
https://github.com/EPFL-LCSB/redhuman/tree/master/redhuman/GEMs/

The Bacillus model was taken from: Henry et al. 2009 Genome Biology , File original named "13059_2009_2219_MOESM4_ESM.xml" (https://link.springer.com/article/10.1186/gb-2009-10-6-r69#MOESM4)
The original supplementary file was splitted into BacillusMetData.mat (original Table S1) and BacillusRxns Data.mat (original Table S2)

The E. Coli model was retrieved: 07/2020 http://bigg.ucsd.edu/models/iJR904
Thermodynamic data was matched to this. Data Taken from: https://www.sciencedirect.com/science/article/pii/S0006349513006851#mmc2 (original Document S2)
For the latter supplementary was reduced to the data needed.

