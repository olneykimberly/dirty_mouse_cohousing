#----------------- Libraries
# lib path .libPaths(c("/R/x86_64-pc-linux-gnu-library", "/opt/R/4.5.2/lib/R/library"),
.libPaths(c("/home/kolney/R/x86_64-pc-linux-gnu-library", "/opt/R/4.5.2/lib/R/library"))
# R version 4.5.2
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(DESeq2) 
library(glmGamPoi)
require(openxlsx)
library(ggrepel)
library(devtools)
library(reshape2)
library(edgeR)  
library(limma) 
library(tidyverse)
library(GenomicFeatures)
library(data.table)
library(gplots)
library(variancePartition)
library(cowplot)
library(UpSetR)
library(ComplexUpset)

#----------------- Create folder outputs
# Create directory structure if it doesn't exist
# This ensures that saveRDS() and saveToPDF() calls don't fail
required_dirs <- c("results_CH/both_sexes/counts",
                   "results_CH/both_sexes/MDS",
                   "results_CH/both_sexes/library",
                   "results_CH/both_sexes/voom",
                   "results_CH/both_sexes/variance",
                   "results_CH/both_sexes/DEGs",
                   "results_CH/both_sexes/PCA",
                   "results_CH/both_sexes/volcano", 
                   "results_CH/both_sexes/gene_expression_checks", 
                   "results_CH/both_sexes/upset",
                
                   # sex_stratified
                   "results_CH/sex_stratified/counts",
                   "results_CH/sex_stratified/MDS",
                   "results_CH/sex_stratified/library",
                   "results_CH/sex_stratified/voom",
                   "results_CH/sex_stratified/variance",
                   "results_CH/sex_stratified/DEGs",
                   "results_CH/sex_stratified/PCA",
                   "results_CH/sex_stratified/volcano", 
                   "results_CH/sex_stratified/gene_expression_checks", 
                   "results_CH/sex_stratified/upset",
                    "results_CH/sex_stratified/log2FC_correlations")

# Loop through and create them
lapply(required_dirs, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))

#----------------- Define variables
typeOfCount <- c("ReadsPerGene.out.tab") 

#----------------- Data
metadata <- read.delim("metadata.tsv", header = TRUE, sep = "\t")
# metadata <- metadata %>% filter(sample != "CFB2") remove outliers here

# Order the groups and samples by group
info_ordered <- metadata %>%
  dplyr::mutate(group = factor(group, levels = c("Clean", "Bedding", "CH"))) %>%
  dplyr::arrange(group)


# Ensure 'sample' is a factor with levels in the desired order (group-wise)
info_ordered$sample <- factor(info_ordered$sample, levels = unique(info_ordered$sample))
metadata <- info_ordered
rm(info_ordered)


#----------------- Functions
saveToPDF <- function(...) {
  d = dev.copy(pdf,...)
  dev.off(d)
}

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

fix_gene_names <- function(dge) {
  
  # replace NA gene_name with gene_id
  dge$genes <- dge$genes %>%
    mutate(gene_name = coalesce(gene_name, gene_id))
  
  # report duplicates BEFORE fixing
  gene_names <- dge$genes$gene_name
  dup_genes <- unique(gene_names[duplicated(gene_names) | duplicated(gene_names, fromLast = TRUE)])
  
  if (length(dup_genes) > 0) {
    message("Duplicated genes found: ", paste(dup_genes, collapse = ", "))
  } else {
    message("No duplicated genes found.")
  }
  
  # make unique
  dge$genes$gene_name <- make.unique(dge$genes$gene_name)
  
  # sanity check
  stopifnot(!any(duplicated(dge$genes$gene_name)))
  
  return(dge)
}

# -----------------------------
# Function to make bidirectional DEG plot
# -----------------------------
make_bidirectional_deg_plot <- function(df, sex_label, qval, lfc.cutoff, outdir) {
  
  # make sure input is a data frame
  df <- as.data.frame(df)
  
  # add missing columns if needed
  if (!"Down" %in% colnames(df)) df$Down <- 0
  if (!"Up" %in% colnames(df)) df$Up <- 0
  
  # keep only Down and Up in the correct order
  df <- df[, c("Down", "Up"), drop = FALSE]
  
  # make rownames exist
  if (is.null(rownames(df))) {
    stop("Input df must have row names corresponding to comparisons.")
  }
  
  df_plot <- df %>%
    rownames_to_column(var = "comparison") %>%
    pivot_longer(
      cols = c("Down", "Up"),
      names_to = "Direction",
      values_to = "n"
    ) %>%
    mutate(
      n = as.numeric(n),
      n_plot = ifelse(Direction == "Down", -n, n),
      comparison = factor(comparison, levels = rev(rownames(df)))
    )
  
  max_n <- max(abs(df_plot$n_plot), na.rm = TRUE)
  x_lim <- max_n * 1.25
  
  p <- ggplot(df_plot, aes(x = comparison, y = n_plot, fill = Direction)) +
    geom_col(width = 0.7) +
    geom_text(
      aes(label = n),
      hjust = ifelse(df_plot$n_plot > 0, -0.15, 1.15),
      size = 3
    ) +
    coord_flip() +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    scale_fill_manual(values = c("Down" = "blue", "Up" = "red")) +
    scale_y_continuous(
      limits = c(-x_lim, x_lim),
      labels = function(x) abs(x),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      title = paste0("Number of DEGs - ", sex_label, "\nq < ", qval, " & |log2FC| > ", lfc.cutoff),
      x = "",
      y = "Number of DEGs"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    )
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  pdf(
    file = paste0(outdir, "/DEGs_bidirectional_", tolower(sex_label), ".pdf"),
    width = 6,
    height = 5
  )
  print(p)
  dev.off()
  
  return(p)
}

pretty_comp <- function(x) {
  gsub("_", " ", x)
}

deg_file <- function(comparison, sex) {
  file.path(deg_dir, paste0(comparison, "_q1.00_", sex, ".txt"))
}

detect_gene_col <- function(df) {
  candidates <- c("gene", "gene_name", "Gene", "SYMBOL", "symbol", "external_gene_name")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

detect_fc_col <- function(df) {
  candidates <- c("log2FoldChange", "avg_log2FC", "avg_logFC", "logFC", "LFC")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_deg_fc <- function(path) {
  if (!file.exists(path)) {
    message("Missing file: ", path)
    return(NULL)
  }

  d <- tryCatch(
    read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )

  if (is.null(d) || nrow(d) == 0) {
    message("Could not read or empty file: ", path)
    return(NULL)
  }

  gene_col <- detect_gene_col(d)
  fc_col   <- detect_fc_col(d)

  if (is.na(gene_col) || is.na(fc_col)) {
    message("Required columns not found in: ", path)
    message("Found columns: ", paste(colnames(d), collapse = ", "))
    return(NULL)
  }

  d <- d[, c(gene_col, fc_col), drop = FALSE]
  colnames(d) <- c("gene", "log2FC")

  d$gene <- as.character(d$gene)
  d$log2FC <- suppressWarnings(as.numeric(d$log2FC))

  d <- d[!is.na(d$gene) & d$gene != "" & !is.na(d$log2FC), , drop = FALSE]

  if (nrow(d) == 0) {
    return(NULL)
  }

  # Deduplicate genes by keeping the largest absolute FC
  if (anyDuplicated(d$gene) > 0) {
    d <- d[order(d$gene, -abs(d$log2FC)), , drop = FALSE]
    d <- d[!duplicated(d$gene), , drop = FALSE]
  }

  return(d)
}

merge_fc <- function(d1, d2, x_name, y_name) {
  if (is.null(d1) || is.null(d2)) return(NULL)

  colnames(d1) <- c("gene", "x")
  colnames(d2) <- c("gene", "y")

  m <- merge(d1, d2, by = "gene", all = FALSE)
  if (nrow(m) == 0) return(NULL)

  m$x_name <- x_name
  m$y_name <- y_name
  m
}

pick_discordant <- function(df, top_n = 5) {
  disc <- df[(df$x > 0 & df$y < 0) | (df$x < 0 & df$y > 0), , drop = FALSE]
  if (nrow(disc) == 0) return(disc)
  disc$score <- abs(df$x[match(disc$gene, df$gene)] - df$y[match(disc$gene, df$gene)])
  disc <- disc[order(-disc$score), , drop = FALSE]
  disc[seq_len(min(top_n, nrow(disc))), , drop = FALSE]
}

pick_concordant <- function(df, top_each = 3) {
  df$score <- abs(df$x) + abs(df$y)

  upup <- df[df$x > 0 & df$y > 0, , drop = FALSE]
  dndn <- df[df$x < 0 & df$y < 0, , drop = FALSE]

  take_top <- function(sub, n) {
    if (nrow(sub) == 0) return(sub)
    sub <- sub[order(-sub$score), , drop = FALSE]
    sub[seq_len(min(n, nrow(sub))), , drop = FALSE]
  }

  rbind(take_top(upup, top_each), take_top(dndn, top_each))
}

make_corr_plot <- function(df,
                           title_line1,
                           title_line2,
                           title_line3,
                           x_lab,
                           y_lab,
                           n_discordant = 5,
                           n_concordant_each = 3) {

  lim <- max(abs(c(df$x, df$y)), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  lim <- lim * 1.05

  r <- suppressWarnings(cor(df$x, df$y, method = "pearson", use = "pairwise.complete.obs"))
  n <- nrow(df)

  title_line4 <- paste0(
    "r = ", ifelse(is.finite(r), sprintf("%.3f", r), "NA"),
    "   n = ", n
  )

  disc_df <- pick_discordant(df, top_n = n_discordant)
  conc_df <- pick_concordant(df, top_each = n_concordant_each)

  if (nrow(disc_df) > 0 && nrow(conc_df) > 0) {
    conc_df <- conc_df[!conc_df$gene %in% disc_df$gene, , drop = FALSE]
  }

  do_repel <- requireNamespace("ggrepel", quietly = TRUE)

  p <- ggplot(df, aes(x = x, y = y)) +
    annotate("rect", xmin = 0, xmax = lim, ymin = 0, ymax = lim,
             fill = "lightpink", alpha = 0.35) +
    annotate("rect", xmin = -lim, xmax = 0, ymin = -lim, ymax = 0,
             fill = "lightblue", alpha = 0.35) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.6, color = "grey40") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50") +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50") +
    geom_point(size = 1.2, color = "black", alpha = 0.8) +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(
      title = paste(title_line1, title_line2, title_line3, title_line4, sep = "\n"),
      x = x_lab,
      y = y_lab
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(size = 9, face = "plain", lineheight = 1.05),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8)
    )

  if (nrow(conc_df) > 0) {
    if (do_repel) {
      p <- p + ggrepel::geom_text_repel(
        data = conc_df,
        aes(label = gene),
        size = 8 / ggplot2::.pt,
        color = "black",
        min.segment.length = 0,
        box.padding = 0.2,
        point.padding = 0.15,
        max.overlaps = 50
      )
    } else {
      p <- p + geom_text(
        data = conc_df,
        aes(label = gene),
        size = 8 / ggplot2::.pt,
        color = "black",
        vjust = -0.6
      )
    }
  }

  if (nrow(disc_df) > 0) {
    if (do_repel) {
      p <- p + ggrepel::geom_text_repel(
        data = disc_df,
        aes(label = gene),
        size = 8 / ggplot2::.pt,
        color = "red3",
        min.segment.length = 0,
        box.padding = 0.2,
        point.padding = 0.15,
        max.overlaps = 50
      )
    } else {
      p <- p + geom_text(
        data = disc_df,
        aes(label = gene),
        size = 8 / ggplot2::.pt,
        color = "red3",
        vjust = -0.6
      )
    }
  }

  p
}

make_missing_plot <- function(title_line1, title_line2, title_line3, msg = "Missing file(s) or no overlapping genes") {
  ggplot() +
    theme_void() +
    ggtitle(paste(title_line1, title_line2, title_line3, msg, sep = "\n")) +
    theme(plot.title = element_text(size = 9, face = "plain", lineheight = 1.05))
}

save_three_panel <- function(p1, p2, p3, out_pdf, width = 13.5, height = 4.2) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    comb <- p1 + p2 + p3 + patchwork::plot_layout(nrow = 1)
    ggsave(out_pdf, comb, width = width, height = height)
    return(invisible(NULL))
  }

  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gr <- gridExtra::arrangeGrob(p1, p2, p3, nrow = 1)
    ggsave(out_pdf, gr, width = width, height = height)
    return(invisible(NULL))
  }

  pdf(out_pdf, width = width, height = height)
  print(p1)
  print(p2)
  print(p3)
  dev.off()
}
