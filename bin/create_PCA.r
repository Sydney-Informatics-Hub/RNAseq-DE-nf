#!/usr/bin/env Rscript

library(dplyr)
library(tidyverse)


# 1. Import data/csv

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

merged_counts_STAR <- args[1]
samples_info <- args[2]

	
countFile <- merged_counts_STAR
raw_cts<-read.table(countFile,header=TRUE,row.names = 1,sep = "\t", check.names=TRUE, stringsAsFactors = FALSE)
head(countFile)

# Actually need to get a transformed version for PCA
countFile <- merged_counts_STAR
trans_cts<-read.table(countFile,header=TRUE,row.names = 1,sep = "\t", check.names=TRUE, stringsAsFactors = FALSE)



designFile <- samples_info
sample_info<-read.table(designFile,header=TRUE,row.names = 1,sep = "\t", stringsAsFactors = FALSE)

#head(sample_info)
sample_info <- cbind(RowNames = rownames(sample_info), sample_info)
rownames(sample_info) <- NULL  # Clear the row names

names(sample_info)[1] <- "sample"

#head(sample_info)


#head (sample_info)

# 2. number of samples is in the "sample_info" table
nrow(sample_info)

# 4. this can be taken from the table of counts
nrow(trans_cts)


#head(raw_cts)

# Exploratory analysis of count data

## Reshape table

# Convert the `raw_cts` table to a "long" format using the `pivot_longer()` function.
# Save it into an object called `raw_cts_long`.


# Convert to long variable
raw_cts_long <- raw_cts %>%
  pivot_longer(everything(), names_to = "sample", values_to = "cts")

# Join with sample information table
raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))

# Make the plot
library(ggplot2)

# head(raw_cts_long)


# Assuming you have already created the ggplot object raw_cts_long


# PCA

## Examine `prcomp()` output
sample_pca <- prcomp(t(trans_cts[, -1]))
pc_scores <- sample_pca$x


## Annotating PC plot

pca_plot <- sample_pca$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample") %>%
  # join with "sample_info" table
  full_join(sample_info, by = "sample") %>%
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(group))) +
  geom_point() +
  # Set labels dynamically based on levels of 'group' variable
  scale_color_discrete(labels = levels(raw_cts_long$group))


# print the result (in this case a ggplot)
# Save the plot to a file
ggsave("pca_plot_STAR.png", pca_plot)




