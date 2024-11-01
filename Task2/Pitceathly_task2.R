
#Script Name: Pitceathly_task2.R
#Author: CMI
#Date: 29/10/2024
#Version: 1.0
#Notes: this is the first task for interview for Pitceathly Lab RA positon
#Purpose: use genes in the Genomics England Panel App and MitoCarta 
        # to create a "super panel" of mitochondrial genes
        # create a BED file with the union of the genes

setwd("/Users/carolinaierardi/Documents/Personal/Currículo/PitceathlyLabRA/Task2")

# Load necessary libraries
library(readr)          # For reading TSV and CSV files
library(readxl)         # For reading Excel files
library(dplyr)          # For data manipulation
library(GenomicRanges)  # For working with genomic ranges
library(tidyverse)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")

library(biomaRt)


# Load the gene lists
panel_genes = read_tsv("Mitochondrial disorders.tsv")
mitocarta_genes = read_excel("Human.MitoCarta3.0.xls", sheet = 2)

print("panel genes have Ensembl ID for both genome builds.
      mitocarta genes have hg19 coordinates")

#Select relevant columns for each panel
#mitocarta only has coordinates for hg19
mitocarta_genes = mitocarta_genes %>%
  dplyr::select(HumanGeneID, Symbol, EnsemblGeneID_mapping_version_20200130,
         hg19_Chromosome, hg19_Start, hg19_Stop)

#Genomics England Panel has Ensembl IDs for both genome builds
panel_genes = panel_genes %>%
  dplyr::select(`Gene Symbol`, HGNC, `EnsemblId(GRch37)`, `EnsemblId(GRch38)`)


#We will use biomaRt library to obtain genome coordinates
mart_hg19 = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="https://grch37.ensembl.org")

obtain_coordinates = function(ensembl_ids, mart) {
  getBM(attributes = c("chromosome_name", 
                       "start_position", 
                       "end_position", 
                       "ensembl_gene_id"), #obtain this from database
        filters = "ensembl_gene_id",       #match based on Ensembl ID
        values = ensembl_ids,              #use list of Ensembl IDs
        mart = mart)                       #use specific genome build
}

format_bed <- function(coordinates) {
  #function to format file for .bed
  coordinates %>%                                                #list of coordinates
    filter(chromosome_name %in% c(1:22, "X", "Y", "MT")) %>%           #only these rows
    mutate(chromosome_name = paste0("chr", chromosome_name)) %>% #format chr column
    #select only the following columns
    dplyr::select(Chr = chromosome_name, 
           Start = start_position, 
           End = end_position, 
           Ensembl_Gene_ID = ensembl_gene_id)
}


#obtain coordinates for each panel
panel_19 = obtain_coordinates(panel_genes$`EnsemblId(GRch37)`, mart_hg19)
mitocarta_19 = obtain_coordinates(mitocarta_genes$EnsemblGeneID_mapping_version_20200130,mart_hg19)

#merge tables for the super panel
merged_genes_19 = panel_19 %>%
  bind_rows(mitocarta_19) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

#format .bed file
hg19_bed = format_bed(merged_genes_19)

#write .bed file
write.table(hg19_bed, "Mitochondrial_Super_Panel_Hg19.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#Now for the newer build
mart_hg38 = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

panel_genes_coordinates_38 = getBM(attributes = c("chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
                                filters = "ensembl_gene_id",
                                values = panel_genes$`EnsemblId(GRch38)`,
                                mart = mart_hg38)

panel_38 = obtain_coordinates(panel_genes$`EnsemblId(GRch38)`,mart_hg38)
mitocarta_38 = obtain_coordinates(mitocarta_genes$EnsemblGeneID_mapping_version_20200130, mart_hg38)


merged_genes_38 = panel_38 %>%
  bind_rows(mitocarta_38) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

hg38_bed = format_bed(merged_genes_38)

write.table(hg38_bed, "Mitochondrial_Super_Panel_Hg38.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


#Visualisations

#Modify data to be in long format with counts for plot
combined_mutations <- hg19_bed %>%
  count(Chr) %>%
  mutate(Source = "hg19") %>%
  bind_rows(
    hg38_bed %>%
      count(Chr) %>%
      mutate(Source = "hg38")
  )


# Plot bar chart
chr_mut = ggplot(combined_mutations, aes(x = Chr, y = n, fill = Source)) +
  geom_bar(stat = "identity", position = "dodge") +
  #change x-axis ticks 
  scale_x_discrete(limits = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
                              "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
                              "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 
                              "chrX", "chrMT"),
                   labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", 
                              "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "MT")) +
  #Add legend
  scale_fill_manual(values = c("hg19" = "darkslategray3", "hg38" = "cornflowerblue")) +
  labs(title = "Mutation per Chromosome",
       x = "Chromosome",
       y = "Count") +
  theme_minimal()

ggsave("chr_mut.png",chr_mut)



ggplot(hg19_bed, aes(x = Chr)) +
  geom_bar(color = "#005f73") +
  scale_x_discrete(limits = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
                              "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
                              "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 
                              "chrX", "chrY", "chrM"),
                   labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", 
                              "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT")) +
  scale_fill_manual(values = c("#0a9396")) +  labs(title = "Mutations in each chromosome",
       x = "Chromosome",
       y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("dna_types.png",dna_types)


