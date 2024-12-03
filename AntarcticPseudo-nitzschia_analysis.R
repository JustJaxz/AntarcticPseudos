title: "Antarctic Pseudo-nitzschia"
author: "JacquiS"
date: "`r Sys.Date()`"  # Automatically inserts the current date

#---- LOAD LIBRARIES ----#
library(devtools)      # Tools for package development
library(phyloseq)      # Analyzing microbiome census data
library(ranacapa)      # For visualizing and analyzing taxonomic data
library(ggplot2)       # Data visualization
library(vegan)         # Community ecology package, useful for diversity analysis
library(dplyr)         # Data manipulation
library(tidyverse)     # A collection of R packages for data science

#---- DEFINE FUNCTIONS ----#
update_taxa_na <- function(taxa_df) {
  # Converts columns to character type
  for (i in 1:9) {
    taxa_df[, i] <- as.character(taxa_df[, i])
  }
  
  taxa_df[is.na(taxa_df)] <- ""  # Replace any missing values with empty strings
  
  # Assigns "Unclassified" labels to missing taxonomic ranks based on available information
  for (i in 1:nrow(taxa_df)) {
    if (taxa_df[i, 2] == "") {
      kingdom <- paste("Unclassified_Domain_", taxa_df[i, 1], sep = "")
      taxa_df[i, 2:9] <- kingdom
    } else if (taxa_df[i, 3] == "") {
      supergroup <- paste("Unclassified_Supergroup_", taxa_df[i, 2], sep = "")
      taxa_df[i, 3:9] <- supergroup
    } else if (taxa_df[i, 4] == "") {
      division <- paste("Unclassified_Division_", taxa_df[i, 3], sep = "")
      taxa_df[i, 4:9] <- division
    } else if (taxa_df[i, 5] == "") {
      subdivision <- paste("Unclassified_SubDivision_", taxa_df[i, 4], sep = "")
      taxa_df[i, 5:9] <- subdivision
    } else if (taxa_df[i, 6] == "") {
      class <- paste("Unclassified_Class_", taxa_df[i, 5], sep = "")
      taxa_df[i, 6:9] <- class
    } else if (taxa_df[i, 7] == "") {
      order <- paste("Unclassified_Order_", taxa_df[i, 6], sep = "")
      taxa_df[i, 7:9] <- order
    } else if (taxa_df[i, 8] == "") {
      family <- paste("Unclassified_Family_", taxa_df[i, 7], sep = "")
      taxa_df[i, 8:9] <- family
    } else if (taxa_df[i, 9] == "") {
      taxa_df$R9[i] <- paste("Unclassified_Genus_", taxa_df$R8[i], sep = "_")
    }
  }
  
  return(taxa_df)  # Returns the updated taxonomic data frame
}

#---- DEFINE ----#
#Sites of interest
polarsites <- c("H1", "H3", "H4", "H5", "H6", "H7", "B01", "A01", "P7", "P8", "P9", "P10", "P11", "P3", "P4", "P5", "P6")
#Taxonomy Divisions for tax_tables
taxo.divisions <- c("Domain", "Supergroup", "Division", "Sub-division", "Class", "Order", "Family", "Genus", "Species")
#Genera of interest
emc.genus <- c("Pseudo-nitzschia")

#-----------------------------------#
#-------- Load & Clean data --------#
#-----------------------------------#
v9 <- readRDS("C:/Users/jacquis/Documents/R/v9.cleaned2.rds")  # Load the cleaned phyloseq object

df <- sample_data(v9)  # Extract sample data

# Clean ASV data
v9 <- filter_taxa(v9, function(x) sum(x) > 0, TRUE)  # Remove taxa with zero counts

# Dataframe
v9.df <- data.frame(tax_table(v9))  # Convert the taxonomic table to a data frame

# Subset the phyloseq object to sites of interest
v9.ice <- subset_samples(v9, siteID %in% polarsites)  # Subset samples to include only polar sites
v9.ice <- filter_taxa(v9.ice, function(x) sum(x) > 0, TRUE)  # Filter out zero-count taxa again

# Calculate total ASVs for each site
total_asvs <- colSums(otu_table(v9.ice))  # Sum the ASVs across each site

# Create subsets of data for specific analysis
v9.ice.emc <- subset_taxa(v9.ice, R8 %in% emc.genus)  # Subset for specific genera (e.g., Pseudo-nitzschia)
v9.ice.emc.tot <- subset_taxa(v9.ice, R5 %in% )  # Subset for specific classes of interest

# Dataframe
v9.ice.df <- data.frame(tax_table(v9.ice.emc))  # Convert the taxonomic table for Pseudo-nitzschia to a data frame

v9.ice.df <- update_taxa_na(v9.ice.df)  # Update missing taxonomic ranks with unclassified labels
colnames(v9.ice.df) <- taxo.divisions  # Rename the columns with the taxonomy divisions
tax_table(v9.ice.emc) <- as.matrix(v9.ice.df)  # Update the tax_table in the phyloseq object
sample_data(v9.ice.emc)  # View the updated sample data


#-----------------------------------#
#--- Calculating ASV proportions ---#
#-----------------------------------#

# Merge samples and calculate ASV counts
v9.ice.pseudo_mrg <- merge_samples(v9.ice.emc, "siteID")  # Merge samples by siteID for Pseudo-nitzschia
v9.ice.emc.tot_mrg <- merge_samples(v9.ice.emc.tot, "siteID")  # Merge samples for the broader EMC classes
v9.ice_mrg <- merge_samples(v9.ice, "siteID")  # Merge all samples by siteID

# Get the sum of ASVs per site
pseudo_asvs <- sample_sums(v9.ice.pseudo_mrg)  # Sum the Pseudo-nitzschia ASVs per site
emc_asvs <- sample_sums(v9.ice.emc.tot_mrg)  # Sum the EMC class ASVs per site
total_asvs <- sample_sums(v9.ice_mrg)  # Sum the total ASVs per site

# Combine into a single data frame
sample_sums_df <- data.frame(
  Site = names(pseudo_asvs),
  Pseudo_ASVs = pseudo_asvs,
  EMC_ASVs = emc_asvs,
  Total_ASVs = total_asvs
)

# Calculate the proportion of Pseudo-nitzschia ASVs to EMC ASVs
sample_sums_df$Pseudo_Proportion <- sample_sums_df$Pseudo_ASVs / sample_sums_df$EMC_ASVs

# Display the data frame
print(sample_sums_df)


#-----------------------------------#
#---------- Plotting data ----------#
#-----------------------------------#

#----------- Bubble plot -----------#
#-- Make new factor to keep year and site separate ---#
# Create a new factor in sample_data combining "siteID" and "year"
sample_data(v9.ice.emc)$site_year <- with(sample_data(v9.ice.emc), paste(siteID, year, sep = "_"))

# Merge samples using the new "site_year" factor
v9.ice.emc_mrg <- merge_samples(v9.ice.emc, "site_year")

# Optionally, view the sample_data to confirm merging
sample_data(v9.ice.emc_mrg)

# Ensure that the sample data includes site information
site_info <- data.frame(sample_data(v9.ice.emc_mrg))

# Create a dataframe that includes taxonomic information for Pseudo-nitzschia
pseudo_nitzschia_taxa <- data.frame(tax_table(v9.ice.emc_mrg))

# Extract the relevant abundance data for Pseudo-nitzschia
pseudo_nitzschia_abundance <- data.frame(otu_table(v9.ice.emc_mrg))

# Combine Abundance with SiteID
pseudo_nitzschia_detection <- cbind(site_info, pseudo_nitzschia_abundance)
pseudo_nitzschia_detection$siteID <- rownames(pseudo_nitzschia_detection)

# Ensure there's a "siteID" column in the site_info
pseudo_nitzschia_detection <- pseudo_nitzschia_detection %>%
  select(siteID, everything())  # Assuming "siteID" is the column name for your sites

# Step 2: Summarize Abundance Data by siteID and ASV
# Melt the dataframe to have siteID, ASV, and Abundance in long format

colnames(pseudo_nitzschia_detection) <- sub(
  pattern = "CGCCCGTCGCACCTACCGATTGAATGGTCCGGTGAAGCCTCGGGATTGTGGCTGGTTTCCTTTATTGGAATCTGACCACGAGAACCTGTCTAAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC", 
  replacement = "ASV1", 
  colnames(pseudo_nitzschia_detection)
)

colnames(pseudo_nitzschia_detection) <- sub(
  pattern = "CGCCCGTCGCACCTACCGATTGAATGGTCCGGTGAAGCCTCGGGATTGTGGTTGGTTTCCTTTATTGGAATCTGACCACGAGAACCTGTCTAAACCTTATCATTTAGAGGAAGGTGAAGTCGTAACAAGGTTTCC", 
  replacement = "ASV2", 
  colnames(pseudo_nitzschia_detection)
)

melted_data <- pseudo_nitzschia_detection %>%
  pivot_longer(
    cols = starts_with("ASV"),   # Assuming your ASV columns start with "ASV"
    names_to = "ASV", 
    values_to = "Abundance"
  )

# Generate the bubble plot
bubble_plot <- ggplot(melted_data, aes(x = siteID, y = ASV, size = Abundance)) +
  geom_point(alpha = 0.6) +
  scale_size_area(max_size = 10) +  # Adjust max_size for bubble scaling
  labs(
    title = "ASV Abundance at Each Site",
    x = "Site ID",
    y = "ASV",
    size = "Abundance"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(bubble_plot)


#------------------------------#
#--- Cape Evans Time Series ---#

# Step 1: Filter the data for sites starting with "H1" and "P5" and for ASV1 only
filtered_data <- melted_data %>%
  filter(grepl("^H1|^P5", siteID) & ASV == "ASV1")

# Step 2: Create a time series plot for ASV1 abundance at each site
timeSeries <- ggplot(filtered_data, aes(x = year, y = Abundance, color = siteID, group = siteID)) +
  geom_line() +           # Line plot for the time series with specified line thickness
  geom_point(size = 3) +          # Points to mark the data points
  labs(
    title = "Pseudo-nitzschia ASV1 Abundance Over Time at Selected Sites",
    x = "Year",
    y = "Abundance",
    color = "Site ID"
  ) +
  theme_minimal() +               # Use a minimal theme as a base
  theme(
    panel.grid = element_blank(),                   # Remove the background grid
    axis.line = element_line(color = "black"),      # Add x and y axis lines
    axis.ticks = element_line(color = "black"),     # Add axis ticks
    axis.ticks.length = unit(0.2, "cm"),            # Set tick marks length
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
  )

print(timeSeries)
  
