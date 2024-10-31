

#Script Name: Pitceathly_task1.R
#Author: CMI
#Date: 26/10/2024
#Version: 1.0
#Notes: this is the first task for interview for Pitceathly Lab RA positon
#Purpose: determine pathogenic mitochondrial variants and
          #create data for clinicians 


#Change environment to import data
setwd("/Users/carolinaierardi/Documents/Personal/Currículo/PitceathlyLabRA/Task1")

library(dplyr)
library(stringr)
library(stringi)
library(gridExtra)
library(ggplot2)


#references
#https://samtools.github.io/hts-specs/VCFv4.2.pdf
#https://stackoverflow.com/questions/32513776/how-to-read-vcf-file-in-r


#Import data
mito_var_names = readLines("task_example.vcf")
mito_var = read.table("task_example.vcf", stringsAsFactors = F)
mitomap = read.csv("ConfirmedMutations  MITOMAP  Foswiki.csv")

mito_var_names = mito_var_names[-(grep("#CHROM",mito_var_names)+1):-(length(mito_var_names))]
vcf_names = unlist(strsplit(mito_var_names[length(mito_var_names)],"\t"))

names(mito_var) = vcf_names

#Data investigation


# extract and plot AF
mito_var <- mito_var %>%
  mutate(AF = as.numeric(str_extract(INFO, "(?<=AF=)[0-9.]+")))

#Make figure of AF distribution

png(filename = "distr_AF.png")
hist(mito_var$AF, main = "Distribution of Heteroplasmy Levels",
     xlab = "AF", col = "seagreen")
dev.off()
unique_samples = length(unique(mito_var$SAMPLE)) # Assuming SampleID is the column name


# Define mutation types by keywords in the "Allele" column
mitomap <- mitomap %>%
  mutate(
    Mutation_Type = case_when(
      str_detect(Allele, "del") ~ "deletion",
      str_detect(Allele, "ins") ~ "insertion",
      str_detect(Allele, ">") ~ "substitution",
      str_detect(Allele, "inv") ~ "inversion",
      TRUE ~ "other"
    ),
    
    # Extract REF and ALT based on mutation type
    REF = case_when(
      Mutation_Type == "substitution" ~ str_extract(Allele, "[ACGT]+(?=>)"),    # Capture sequence before >
      Mutation_Type == "deletion" ~ str_extract(Allele, "(?<=del)[ACGT]+"),     # Capture sequence after 'del'
      Mutation_Type == "insertion" ~ "",                                        # For insertions, REF is empty
      Mutation_Type == "inversion" ~ str_extract(Allele, "[ACGT]+(?=inv)"),    # Capture sequence after 'inv'
      TRUE ~ NA_character_
    ),
    
    ALT = case_when(
      Mutation_Type == "substitution" ~ str_extract(Allele, "(?<=>)[ACGT]+"),   # Capture sequence after >
      Mutation_Type == "insertion" ~ str_extract(Allele, "(?<=ins)[ACGT]+"),    # Capture sequence after 'ins'
      Mutation_Type == "deletion" ~ "*",                                         # For deletions, ALT is empty
      Mutation_Type == "inversion" ~ stri_reverse(REF),                         # Reverse the REF sequence for ALT in inversions
      TRUE ~ NA_character_
    )
  )


#find the variants that are in both tables
pathogenic_variants <- mito_var %>%
  inner_join(mitomap, by = c("POS" = "Position", "REF", "ALT"))

#make a clinician table with the information requested
clinician_table <- pathogenic_variants %>%
  select(SAMPLE, POS, REF, ALT, AF) %>%
  arrange(POS)

pdf("clinician_table.pdf")
grid.table(clinician_table)
dev.off()

#Visualisations
                ## Make plot of types of mutation
# Categorize mutation types in VCF file
mito_var <- mito_var %>%
  mutate(Mutation_Type = case_when(
    nchar(REF) == 1 & nchar(ALT) == 1 ~ "substitution",
    nchar(REF) < nchar(ALT) ~ "insertion",
    nchar(REF) > nchar(ALT) ~ "deletion",
    TRUE ~ "Other"
  ))

# Combine mutation types from both files
combined_mutation_types <- mito_var %>%
  count(Mutation_Type) %>%
  mutate(Source = "VCF") %>%
  bind_rows(
    mitomap %>%
      count(Mutation_Type) %>%
      mutate(Source = "MitoMap")
  )


# Plot bar chart
mut_types = ggplot(combined_mutation_types, aes(x = Mutation_Type, y = n, fill = Source)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("VCF" = "#005f73", "MitoMap" = "#94d2bd")) +
  labs(title = "Types of Mutation",
       x = "Mutation Type",
       y = "Count") +
  theme_minimal()
ggsave("mut_types.png",mut_types)

                  ##Make plot of mutation position

vcf_pos = ggplot(mito_var, aes(x = POS)) +
  geom_histogram(binwidth = 100, fill = "#94d2bd", color = "#005f73", alpha = 0.7) +
  labs(title = "Mutation Positions in VCF Data",
       x = "Position in Mitochondrial Genome",
       y = "Count") +
  theme_minimal()
ggsave("vcf_pos.png",vcf_pos)

# Plot histogram for MitoMap data positions
mito_pos = ggplot(mitomap, aes(x = Position)) +
  geom_histogram(binwidth = 100, fill = "#0a9396", color = "#005f73", alpha = 0.7) +
  labs(title = "Mutation Positions in MitoMap Data",
       x = "Position in Mitochondrial Genome",
       y = "Count") +
  theme_minimal()
ggsave("mito_pos.png",mito_pos)


              ## Bar chart of DNA types in MitoMap
dna_types = ggplot(mitomap, aes(x = Locus.Type, fill = Locus.Type)) +
  geom_bar(color = "#005f73") +
  scale_fill_manual(values = c("#0a9396", "#005f73", "#94d2bd", "#ee9b00", "#ca6702")) +
  labs(title = "DNA Types Affected by Mutations (MitoMap)",
       x = "DNA Type",
       y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("dna_types.png",dna_types)


## Position vs Pathgenicity

mito_var$Pathgenic = ifelse(mito_var$POS %in% pathogenic_variants$POS, "Pathogenic", "Non-Pathogenic")

# Plot the scatter plot of Position vs. Heteroplasmy Level
position_patho = ggplot(mito_var, aes(x = POS, y = AF, color = Pathgenic)) +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = c("Pathogenic" = "#ee9b00", "Non-Pathogenic" = "#0a9396")) +
  labs(title = "Pathogenicity vs. Position",
       x = "Position",
       y = "Heteroplasmy Level",
       color = "Pathogenicity") +
  theme_minimal()

ggsave("position_patho.png",position_patho)
