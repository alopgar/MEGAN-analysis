# MEGAN-analysis
DIAMOND + MEGAN pipeline for Nanopore MinION metagenome data.

## 0. SOFTWARE REQUIREMENTS:
For the correct functioning of this scripts, the installation of several software is required:
..- **prinseq-lite**: https://sourceforge.net/projects/prinseq/files/standalone/
..- **DIAMOND**: Currently available at CESGA server modules
..- **MEGAN6**: http://ab.inf.uni-tuebingen.de/software/megan6/
..- **SPARCC**: Inside a .conda environment. https://bitbucket.org/yonatanf/sparcc/src/default/
..- **R**: With packages tidyr and PCIT

## 1. Script execution:
a) Use **sh MinION_MEGAN.sh -h**  
b) Change variable paths inside *MinION_MEGAN.sh*  
c) Use **sbatch MinION_MEGAN.sh -i** for running the script.  

## 2. Getting trimming and performance statistics:
Use **sbatch MinION_MEGAN.sh -F**  
Output: All files will be stored in output directory ($WDIR)  
- Minion_trim.tsv -> Trimming summary
- Minion_dmnd.tsv -> Running times for DIAMOND
- Minion_megan.tsv -> Running times for MEGAN

## 3. Final processing:
We recommend execunting in console (**sh MinION_MEGAN_finalprocessing.sh**), as this is an R script which might be modified as we run it. 
- **_First step_**: Create final taxonomy file using all samples MEGAN files. It is needed to check which phyla are from plant or animal origin and remove them.
- **_Second step_**: create input file for sparcc.

## 4. SPARCC and networks files:
- **_Input_**: .txt (tsv) table formatted in **R** from *rma_taxa_count_filtered.csv*
- **_Output_**: Correlation matrix

After, R processing of the correlation matrix is performed in order to get an edgelist for cytoscape or other network software.
