# Load packages
library(tidyverse)
library(janitor)

# Step 1: Data
#   - Here we load the data from PEAKS and perform some standardization steps as contents may differ slightly file-to-file.
#   - This code is written to operate on the "protein-peptides.csv" file that is exported from PEAKS 10pro after a database search. 

source_file <- read.csv("Path to your file here") %>% # update file path as needed
  clean_names()

#   - This line renames the "Area" column as the generalized "intensity" for simplification. 
#     Check your source file to ensure the area data is in column 13 before beginning. 
colnames(source_file)[13] <- "sample_intensity"


source_file <- source_file %>%
  mutate_at(vars(protein_accession, peptide), as.character) %>%
  mutate_at(vars(sample_intensity, start, end), as.numeric) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(across(where(is.character), ~ ifelse(. %in% c("", " "), NA, .))) %>%
  mutate(across(where(is.character), tolower)) %>% 
  as_tibble()

#   -This line extracts the name of the PEAKS source file from the .csv file and includes it as an identifier in your results.
identifier <- source_file %>% slice(1) %>% pull(source_file)

glimpse(source_file)


#   - For both serum albumin and IgG, the accession numbers are based upon the Uniprot Consortium database. 
#   - Start and end values are determined by the location of the most abundant relevant peptide. 


# Step 2: Serum Albumin
#   - These steps filter your source file for the protein of interest, 
#     locates peptides occurring at the location of the most prevalent amino acid sequence (pre-determined), 
#     selects those with deamidated regions, and performs the ratio calculation. 

serum_albumin <- source_file %>% 
  filter(protein_accession == "p02768") %>%
  filter(start >= 438, end <= 456) %>% 
  filter(str_detect(peptide, "pq") | str_detect(peptide, "rn"))

serum_albumin_sum_total <- sum(serum_albumin$sample_intensity)
serum_albumin_E <- serum_albumin %>% filter(str_detect(peptide, "98")) 
serum_albumin_E_sum_total <- sum(serum_albumin_E$sample_intensity)

# Calculate SA deamidation ratio
serum_albumin_percent_deamidation <- (serum_albumin_E_sum_total/serum_albumin_sum_total)*100



# Step 3: IgG1
#   - Same process as above, for a second protein of interest.

igg1 <- source_file %>%
  filter(protein_accession == "p01857") %>%
  filter(start >= 184, end <= 201) %>%
  filter(str_detect(peptide, "hq") | str_detect(peptide, "ln"))

igg1_sum_total <- sum(igg1$sample_intensity)
igg1_E <- igg1 %>% filter(str_detect(peptide, "98"))
igg1_E_sum_total <- sum(igg1_E$sample_intensity)

# Calculate IgG deamidation ratio
igg1_percent_deamidation <- (igg1_E_sum_total/igg1_sum_total)*100



# Step 4: view results

results <- tibble(`Serum Albumin` = serum_albumin_percent_deamidation,
                  `IgG1` = igg1_percent_deamidation,
                  `Identifier` = identifier)
results




