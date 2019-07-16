#!/usr/bin/env Rscript

####################################################
## FINAL TAXONOMY/RA FILE (IN R)
####################################################
module load cesga/2018 gcc/6.4.0 R/3.5.1
R
library(tidyr)

WD <- "/mnt/lustre/scratch/home/otras/ini/alg/Results/Metagenome/METALGEN/"
infofiles <- list.files(paste0(WD, "4_files.rmainfo/"), pattern="rma.RA.txt")
infofiles.dir <- paste0(WD, "4_files.rmainfo/", infofiles)

infolist <- lapply(infofiles.dir, read.table, sep="\t", stringsAsFactors = F)
names(infolist) <- gsub("_trim.*$", "", infofiles)
for (i in seq(infolist)){names(infolist[[i]]) <- c("Rank", "TaxPath", names(infolist[i]))}

rma_merged <- Reduce(function(x, y) merge(x, y, all=TRUE), infolist)
rma_merged[is.na(rma_merged)] <- 0
rma_merged <- separate(data = rma_merged, col = TaxPath, sep = "; ",
                       into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
rma_merged[] <- lapply(rma_merged, gsub, pattern=';', replacement='')

write.csv(rma_merged, file = paste0(WD,"4_files.rmainfo/","rma_taxa_count.csv"), quote = F, row.names = F)

euk <- rma_merged[grepl("Eukaryota", rma_merged$Kingdom),]
euk <- euk[order(euk$Phylum),]
dput(gsub(unique(euk$Phylum), pattern = "\\[P\\] ", replacement = ""))

#Create vector to remove metazoa and plantae. Then filter:
nope <- c("Annelida", "Arthropoda", "Chordata", "Cnidaria", "Mollusca", "Nematoda", "Streptophyta")
rma_filtered <- rma_merged[!grepl(paste(nope, collapse="|"), rma_merged$Phylum),]

#Check present Eukaryota (after filtering by 'nope'):
rma_euk_filt <- rma_filtered[grepl("Eukaryota", rma_filtered$Kingdom),]

rma_filtered[is.na(rma_filtered)] <- ""
rma_filtered[rma_filtered$Kingdom=="","Kingdom"] <- "[SK] unknown"
write.csv(rma_filtered[,-1], file = paste0(WD,"4_files.rmainfo/","rma_taxa_count_filtered.csv"), quote = F, row.names = F)

# SPARCC input table:

fill_phy_tax <- function(data, physeq=TRUE){
  prefix <- c("[SK]", "[P]", "[C]", "[O]", "[F]", "[G]", "[S]")
  patt <- "\\[.+\\]\\s"
  if(isTRUE(physeq)){tax.clean <- data.frame(tax_table(data))}
  else {tax.clean <- data}
  for (i in 1:7){tax.clean[,i] <- as.character(tax.clean[,i])}
  tax.clean[is.na(tax.clean)] <- ""
  for (i in 1:nrow(tax.clean)){
    if (tax.clean[i,2] == ""){
      kingdom <- paste0("Kingdom_", gsub(pattern = patt, replacement = "", tax.clean[i,1]))
      for (j in 2:7){
        tax.clean[i, j] <- paste0(prefix[1], " ", kingdom)
      }
    }
	else if (tax.clean[i,3] == ""){
      phylum <- paste0("Phylum_", gsub(pattern = patt, replacement = "", tax.clean[i,2]))
      for (j in 3:7){
        tax.clean[i, j] <- paste0(prefix[2], " ", phylum)
      }
    }
	else if (tax.clean[i,4] == ""){
      class <- paste0("Class_", gsub(pattern = patt, replacement = "", tax.clean[i,3]))
      for (j in 4:7){
        tax.clean[i, j] <- paste0(prefix[3], " ", class)
      }
    }
	else if (tax.clean[i,5] == ""){
      order <- paste0("Order_", gsub(pattern = patt, replacement = "", tax.clean[i,4]))
      for (j in 5:7){
        tax.clean[i, j] <- paste0(prefix[4], " ", order)
      }
    }
	else if (stringr::str_detect(tax.clean[i,5], "uncultured")){
      order <- paste0("Order_", gsub(pattern = patt, replacement = "", tax.clean[i,4]), "_uncult")
      for (j in 5:7){
        tax.clean[i, j] <- paste0(prefix[5], " ", order)
      }
    }
	else if (tax.clean[i,6] == ""){
      family <- paste0("Family_", gsub(pattern = patt, replacement = "", tax.clean[i,5]))
      for (j in 6:7){
        tax.clean[i, j] <- paste0(prefix[6], " ", family)
      }
    }
	else if (stringr::str_detect(tax.clean[i,6], "uncultured")){
      family <- paste0("Family_", gsub(pattern = patt, replacement = "", tax.clean[i,5]), "_uncult")
      for (j in 6:7){
        tax.clean[i, j] <- paste0(prefix[j], " ", family)
      }
    }
	else if (tax.clean[i,7] == ""){
      tax.clean[i,7] <- paste("[G] Genus", gsub(pattern = patt, replacement = "", tax.clean[i,6]), 
                              sep = "_")
    }
  }
  if(isTRUE(physeq)){tax.return <- tax_table(tax.clean)}
  else {tax.return <- tax.clean}
  return(tax.return)
}


rma_imported <- read.csv(paste0(WD,"4_files.rmainfo/","rma_taxa_count_filtered.csv"), header = T, stringsAsFactors = F, na.strings = c("", "NA"))

rma_sparcc <- fill_phy_tax(rma_imported, FALSE)
rma_sparcc <- rma_sparcc[,-c(1:6)]
write_table(rma_sparcc, file = paste0(WD,"7_networks/","sparcctable.txt"), quote = F, sep = '\t', row.names = F)

quit()
n
module purge
