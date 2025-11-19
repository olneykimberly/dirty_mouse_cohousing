#----------------- Libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(DESeq2) 
require(openxlsx)
library(ggrepel)
library(glmGamPoi)
library(devtools)
library(reshape2)
library(edgeR)  
library(limma)  
library(tidyverse)
library(GenomicFeatures)
library(data.table)
library(philentropy)
library(gplots)
library(variancePartition)
library(NatParksPalettes) # colors
library(cowplot)

#----------------- Define variables
typeOfCount <- c("ReadsPerGene.out.tab") 

#----------------- Functions
saveToPDF <- function(...) {
  d = dev.copy(pdf,...)
  dev.off(d)
}

### --- Function to Extract Gene Expression ---
# Helper function to extract CPM for a specific gene and merge it into the info dataframe
extract_gene_cpm <- function(gene_name, cpm_data, info_data) {
  # Subset for the specific gene
  gene_cpm <- cpm_data %>%
    dplyr::filter(gene_name == !!gene_name)
  
  # Melt and clean the data
  gene_cpm_melt <- gene_cpm %>%
    reshape2::melt() %>%
    dplyr::rename(sample = variable) %>%
    dplyr::filter(sample != "gene_name") %>% # Remove the gene_name column row
    dplyr::select(sample, value) %>%
    dplyr::rename(!!gene_name := value)
  
  # Merge the new gene column into the info dataframe
  info_data <- info_data %>%
    dplyr::left_join(gene_cpm_melt, by = "sample")
  
  return(info_data)
}

fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

makePaddedDataFrame <- function(l, ...) {
  maxlen <- max(sapply(l, length))
  data.frame(lapply(l, na.pad, len = maxlen), ...)
}
#----------------- Data
metadata <- read.delim("/tgen_labs/jfryer/kolney/dirty_mice/dirty_mouse_cohousing/metadata.tsv", header = TRUE, sep = "\t")
# Order the groups and samples by group
info_ordered <- metadata %>%
  dplyr::mutate(group = factor(group, levels = c("Clean", "Bedding", "CH"))) %>%
  dplyr::arrange(group)

# Ensure 'sample' is a factor with levels in the desired order (group-wise)
info_ordered$sample <- factor(info_ordered$sample, levels = unique(info_ordered$sample))
metadata <- info_ordered
rm(info_ordered)
