library(shiny)
library(ggplot2)
library(dplyr)
library(DT)
library(plotly)
library(ggrepel)
library(stats)
library(tidyr)
library(stringr)

# ---------------------------------------------------------
# 1. GLOBAL CONFIGURATION AND DATA LOADING
# ---------------------------------------------------------

# setwd("/tgen_labs/jfryer/kolney/dirty_mice/dirty_mouse_cohousing/scripts/shiny")

projectID <- "CH"
comparisons_list <- c("CH_vs_Bedding", "Bedding_vs_Clean", "CH_vs_Clean")

# Color scales for DT tables
brks <- seq(-1.6, 1.6, .2)
clrs <- colorRampPalette(c("#6baed6", "white", "red"))(length(brks) + 1)

brks2 <- seq(0, 0.05, 0.01)
clrs2 <- colorRampPalette(c("green", "lightgreen", "white"))(length(brks2) + 1)

# ---------------------------------------------------------
# Read metadata
# ---------------------------------------------------------
metadata <- read.delim("metadata.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Standardize sex column
if ("sex" %in% colnames(metadata)) {
  metadata$sex <- toupper(trimws(as.character(metadata$sex)))
} else if ("Sex" %in% colnames(metadata)) {
  metadata$sex <- toupper(trimws(as.character(metadata$Sex)))
} else {
  stop("metadata.tsv must contain a column named 'sex' or 'Sex'.")
}

metadata$sex <- dplyr::case_when(
  metadata$sex %in% c("M", "MALE") ~ "male",
  metadata$sex %in% c("F", "FEMALE") ~ "female",
  TRUE ~ NA_character_
)

if (any(is.na(metadata$sex))) {
  warning("Some values in metadata$sex could not be mapped to 'male' or 'female'.")
}

# Standardize group ordering
metadata <- metadata %>%
  mutate(
    group = factor(group, levels = c("Clean", "Bedding", "CH"))
  ) %>%
  arrange(group)

metadata$sample <- factor(metadata$sample, levels = unique(metadata$sample))

# ---------------------------------------------------------
# Read DEG tables
# ---------------------------------------------------------

# Combined/all-samples DEG results
DEGs_all_raw <- read.delim("DEGs_AllComparisons.txt", stringsAsFactors = FALSE)

DEGs_all <- DEGs_all_raw %>%
  mutate(
    comparison = gsub(" ", "_", comparison),
    comparison = as.character(comparison),
    gene_name = as.character(gene_name),
    gene_id = as.character(gene_id),
    logFC = as.numeric(logFC),
    P.Value = as.numeric(P.Value),
    adj.P.Val = as.numeric(adj.P.Val),
    AveExpr = as.numeric(AveExpr),
    sex = "both"
  ) %>%
  filter(!is.na(logFC), !is.na(P.Value), !is.na(adj.P.Val))

# Sex-stratified DEG results
DEGs_sex_raw <- read.delim("DEGs_AllComparisons_sex_stratified.txt", stringsAsFactors = FALSE)

if (!"sex" %in% colnames(DEGs_sex_raw)) {
  stop("DEGs_AllComparisons_sex_stratified.txt must contain a column named 'sex'.")
}

DEGs_sex <- DEGs_sex_raw %>%
  mutate(
    comparison = gsub(" ", "_", comparison),
    comparison = as.character(comparison),
    sex = case_when(
  toupper(trimws(as.character(sex))) %in% c("M", "MALE") ~ "male",
  toupper(trimws(as.character(sex))) %in% c("F", "FEMALE") ~ "female",
  TRUE ~ tolower(trimws(as.character(sex)))
),
    gene_name = as.character(gene_name),
    gene_id = as.character(gene_id),
    logFC = as.numeric(logFC),
    P.Value = as.numeric(P.Value),
    adj.P.Val = as.numeric(adj.P.Val),
    AveExpr = as.numeric(AveExpr)
  ) %>%
  filter(!is.na(logFC), !is.na(P.Value), !is.na(adj.P.Val))

# For DEG table tab, keep one combined table that includes both + sex-stratified
DEGs_table <- bind_rows(
  DEGs_all %>% mutate(analysis = "both"),
  DEGs_sex %>% mutate(analysis = "sex_stratified")
)

# ---------------------------------------------------------
# Read CPM tables
# ---------------------------------------------------------
cpm_full_all <- read.delim("filtered_cpm_gene_names.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cpm_full_male <- read.delim("filtered_cpm_male_gene_names.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cpm_full_female <- read.delim("filtered_cpm_female_gene_names.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Gene options
gene_options <- sort(unique(c(
  DEGs_all$gene_name,
  DEGs_sex$gene_name
)))
gene_options <- gene_options[!is.na(gene_options) & gene_options != ""]

# ---------------------------------------------------------
# 2. HELPER FUNCTIONS
# ---------------------------------------------------------

get_deg_dataset <- function(sex_choice, DEGs_all, DEGs_sex) {
  if (sex_choice == "both") {
    return(DEGs_all)
  } else {
    return(DEGs_sex %>% filter(sex == sex_choice))
  }
}

safe_hadjpval <- function(data, qval_fixed = 0.05) {
  sig_pvals <- data$P.Value[data$adj.P.Val < qval_fixed & !is.na(data$P.Value)]
  sig_pvals <- sig_pvals[sig_pvals > 0]
  
  if (length(sig_pvals) == 0) {
    return(-log10(qval_fixed))
  } else {
    return(-log10(max(sig_pvals, na.rm = TRUE)))
  }
}

load_gene_counts <- function(gene_name, metadata, cpm_full, sex_subset = NULL) {
  
  if (!("gene_name" %in% colnames(cpm_full))) {
    return(NULL)
  }
  
  if (!(gene_name %in% cpm_full$gene_name)) {
    return(NULL)
  }
  
  gene_row <- cpm_full %>% filter(gene_name == !!gene_name)
  
  exclude_cols <- c("gene_id", "Chr", "width", "gene_name", "gene_type", "type", "GENEID")
  exclude_cols <- intersect(exclude_cols, colnames(gene_row))
  
  sample_cpm <- gene_row %>%
    select(-all_of(exclude_cols)) %>%
    pivot_longer(cols = everything(), names_to = "sample", values_to = "cpm")
  
  data_cpm <- sample_cpm %>%
    left_join(metadata, by = "sample")
  
  data_cpm$gene_name <- gene_name
  data_cpm$cpm <- as.numeric(data_cpm$cpm)
  
  if (!is.null(sex_subset)) {
    data_cpm <- data_cpm %>% filter(tolower(sex) == tolower(sex_subset))
  }
  
  data_cpm <- data_cpm %>%
    filter(!is.na(group), !is.na(cpm))
  
  if (nrow(data_cpm) == 0) {
    return(NULL)
  }
  
  data_cpm$group <- factor(data_cpm$group, levels = c("Clean", "Bedding", "CH"))
  data_cpm <- data_cpm %>% arrange(group)
  data_cpm$sample <- factor(data_cpm$sample, levels = unique(as.character(data_cpm$sample)))
  
  return(data_cpm %>% select(sample, cpm, gene_name, group, sex))
}

make_counts_subplot <- function(data_cpm, gene_name, title_prefix = "All samples") {
  
  if (is.null(data_cpm) || nrow(data_cpm) == 0) {
    return(
      plot_ly() %>%
        layout(title = paste(title_prefix, "-", gene_name, ": no cpm data available"))
    )
  }
  
  custom_colors <- c("Clean" = "skyblue", "Bedding" = "orange", "CH" = "brown")
  
  p_box <- plot_ly(
    data_cpm,
    y = ~cpm,
    x = ~group,
    text = ~paste(
      "Sample:", sample,
      "<br>Group:", group,
      "<br>Sex:", sex,
      "<br>cpm:", round(cpm, 2)
    ),
    color = ~group,
    colors = custom_colors,
    type = "box",
    boxpoints = "all",
    jitter = 0.2,
    marker = list(size = 8),
    showlegend = FALSE
  ) %>%
    layout(
      title = paste(title_prefix, "(Box Plot)"),
      xaxis = list(title = "", showticklabels = FALSE),
      yaxis = list(title = "Counts Per Million (cpm)", zeroline = FALSE)
    )
  
  p_bar <- plot_ly(
    data_cpm,
    x = ~sample,
    y = ~cpm,
    type = "bar",
    color = ~group,
    colors = custom_colors,
    text = ~paste(
      "Sample:", sample,
      "<br>Group:", group,
      "<br>Sex:", sex,
      "<br>cpm:", round(cpm, 2)
    ),
    hoverinfo = "text",
    showlegend = TRUE
  ) %>%
    layout(
      title = paste(title_prefix, "(Individual Sample Values)"),
      xaxis = list(title = "Sample", tickangle = 45),
      yaxis = list(title = "Counts Per Million (cpm)", zeroline = TRUE),
      margin = list(b = 100)
    )
  
  subplot(
    p_box, p_bar,
    nrows = 2,
    heights = c(0.4, 0.6),
    shareX = FALSE,
    titleY = TRUE
  ) %>%
    layout(
      title = list(text = paste(title_prefix, "- Gene:", gene_name)),
      hovermode = "closest"
    )
}

# ---------------------------------------------------------
# 3. UI
# ---------------------------------------------------------

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .container-fluid { padding: 30px; }
      .well { background-color: #f7f7f7; border: 1px solid #e3e3e3; border-radius: 8px; padding: 15px; }
      .sidebar-panel { background-color: #ffffff; padding: 20px; border-right: 1px solid #eee; }
      .main-panel { padding: 20px; }
    "))
  ),
  
  titlePanel("DEG analysis among clean, bedding, and cohoused (CH) mice"),
  h5("Data generated by Donna Roscoe & shiny app created by Kimberly Olney, PhD"),
  h5("Lab of Dr. John Fryer at TGen Arizona"),
  
  tabsetPanel(
    id = "main_tabs",
    
    # -----------------------------------------------------
    # Volcano Tab
    # -----------------------------------------------------
    tabPanel(
      "Volcano Plot",
      sidebarLayout(
        sidebarPanel(
          class = "sidebar-panel",
          width = 3,
          
          selectInput(
            "volcano_comparison",
            "Select Comparison:",
            choices = comparisons_list
          ),
          
          selectInput(
            "volcano_sex",
            "Select Sex:",
            choices = c("Both sexes" = "both", "male" = "male", "female" = "female"),
            selected = "both"
          ),
          
          uiOutput("gene_name_selector_volcano"),
          
          hr(),
          tags$b("Significance threshold:"),
          tags$ul(
            tags$li("q-value < 0.05"),
            tags$li("|Log2FC| > 0.0")
          )
        ),
        
        mainPanel(
          class = "main-panel",
          width = 9,
          h3(textOutput("volcano_title")),
          plotlyOutput(outputId = "volcano", height = "650px")
        )
      )
    ),
    
    # -----------------------------------------------------
    # DEG Table Tab
    # -----------------------------------------------------
    tabPanel(
      "DEG Tables",
      div(
        class = "well",
        fluidRow(
          column(
            width = 3,
            selectInput(
              "comparison_filter",
              "Comparison:",
              c("All", unique(as.character(DEGs_table$comparison)))
            )
          ),
          column(
            width = 3,
            selectInput(
              "table_sex_filter",
              "Sex:",
              c("All", "both", "male", "female")
            )
          ),
          column(
            width = 3,
            selectInput(
              "gene_name_filter",
              "Gene Name:",
              c("All", gene_options)
            )
          ),
          column(
            width = 3,
            selectInput(
              "analysis_filter",
              "Analysis:",
              c("All", unique(DEGs_table$analysis))
            )
          )
        ),
        h4("All Differential Expression Results"),
        DT::dataTableOutput("table")
      ),
      style = "padding: 20px;"
    ),
    
    # -----------------------------------------------------
    # CPM Counts Tab
    # -----------------------------------------------------
    tabPanel(
      "cpm Counts",
      div(
        class = "well",
        fluidRow(
          column(
            width = 3,
            h4("Counts Per Million (cpm) Plots"),
            selectInput(
              "counts_labels",
              "Gene Name:",
              choices = gene_options,
              selected = if ("Lcn2" %in% gene_options) "Lcn2" else gene_options[1]
            )
          ),
          column(
            width = 9,
            h4("All samples"),
            plotlyOutput(outputId = "counts_overall", height = "650px"),
            br(),
            h4("Male"),
            plotlyOutput(outputId = "counts_male", height = "650px"),
            br(),
            h4("Female"),
            plotlyOutput(outputId = "counts_female", height = "650px")
          )
        )
      )
    )
  )
)

# ---------------------------------------------------------
# 4. SERVER
# ---------------------------------------------------------

server <- function(input, output, session) {
  
  data_objects <- reactive({
    list(
      DEGs_all = DEGs_all,
      DEGs_sex = DEGs_sex,
      DEGs_table = DEGs_table,
      metadata = metadata,
      cpm_full_all = cpm_full_all,
      cpm_full_male = cpm_full_male,
      cpm_full_female = cpm_full_female
    )
  })
  
  # -----------------------------------------------------
  # Volcano reactive data
  # -----------------------------------------------------
  volcano_data <- reactive({
    req(input$volcano_comparison, input$volcano_sex)
    
    dat <- get_deg_dataset(
      sex_choice = input$volcano_sex,
      DEGs_all = data_objects()$DEGs_all,
      DEGs_sex = data_objects()$DEGs_sex
    )
    
    dat <- dat %>%
      filter(comparison == input$volcano_comparison)
    
    req(nrow(dat) > 0)
    
    qval_fixed <- 0.05
    lfc_cutoff_fixed <- 0.0
    
    dat <- dat %>%
      mutate(
        Expression = case_when(
          adj.P.Val < qval_fixed & logFC > lfc_cutoff_fixed ~ "Upregulated",
          adj.P.Val < qval_fixed & logFC < -lfc_cutoff_fixed ~ "Downregulated",
          TRUE ~ "Not Significant"
        )
      )
    
    counts <- dat %>%
      group_by(Expression) %>%
      summarise(count = n(), .groups = "drop")
    
    upregulated_count <- counts$count[counts$Expression == "Upregulated"]
    downregulated_count <- counts$count[counts$Expression == "Downregulated"]
    
    list(
      data = dat,
      upregulated = ifelse(length(upregulated_count) > 0, upregulated_count, 0),
      downregulated = ifelse(length(downregulated_count) > 0, downregulated_count, 0)
    )
  })
  
  output$gene_name_selector_volcano <- renderUI({
    selectizeInput(
      "selected_gene_name_volcano",
      "Search/Highlight Gene:",
      choices = c("None" = "None", gene_options),
      selected = if ("Lcn2" %in% gene_options) "Lcn2" else "None",
      options = list(maxOptions = 10, placeholder = "Search or select a gene")
    )
  })
  
  output$volcano_title <- renderText({
    sex_label <- switch(
      input$volcano_sex,
      "both" = "Both sexes",
      "male" = "Male",
      "female" = "Female"
    )
    
    paste(
      "Volcano Plot for:",
      gsub("_vs_", " vs ", input$volcano_comparison),
      "-",
      sex_label
    )
  })
  
  output$volcano <- renderPlotly({
    volcano_list <- volcano_data()
    data <- volcano_list$data
    up_count <- volcano_list$upregulated
    down_count <- volcano_list$downregulated
    
    req(nrow(data) > 0)
    
    qval_fixed <- 0.05
    lfc_cutoff_fixed <- 0.0
    hadjpval <- safe_hadjpval(data, qval_fixed = qval_fixed)
    
    data <- data %>%
      mutate(
        plot_text = paste(
          "Gene:", gene_name,
          "<br>Log2FC:", round(logFC, 3),
          "<br>P.Value:", format.pval(P.Value, digits = 4),
          "<br>Adj.P.Val:", format.pval(adj.P.Val, digits = 4),
          "<br>Comparison:", comparison,
          "<br>Sex:", ifelse("sex" %in% colnames(data), sex, "both"),
          "<br>Status:", Expression
        ),
        is_selected = gene_name == input$selected_gene_name_volcano
      )
    
    max_y <- max(-log10(pmax(data$P.Value, .Machine$double.xmin)), na.rm = TRUE)
    max_x <- max(data$logFC, na.rm = TRUE)
    min_x <- min(data$logFC, na.rm = TRUE)
    
    p <- ggplot(data, aes(x = logFC, y = -log10(P.Value), text = plot_text)) +
      geom_point(
        data = filter(data, !is_selected),
        aes(color = Expression),
        size = 1,
        alpha = 0.6
      ) +
      geom_hline(yintercept = hadjpval, colour = "black", linetype = "dashed") +
      geom_vline(
        xintercept = c(lfc_cutoff_fixed, -lfc_cutoff_fixed),
        colour = "black",
        linetype = "dashed"
      ) +
      geom_text(
        x = max_x,
        y = max_y * 0.95,
        label = paste("Up:", up_count),
        color = "red",
        size = 5,
        hjust = 1,
        vjust = 1
      ) +
      geom_text(
        x = min_x,
        y = max_y * 0.95,
        label = paste("Down:", down_count),
        color = "blue",
        size = 5,
        hjust = 0,
        vjust = 1
      ) +
      {
        if (input$selected_gene_name_volcano != "None" && any(data$is_selected)) {
          geom_point(
            data = filter(data, is_selected),
            color = "black",
            size = 4,
            shape = 21,
            stroke = 1.5,
            fill = "gold"
          )
        }
      } +
      {
        if (input$selected_gene_name_volcano != "None" && any(data$is_selected)) {
          geom_text_repel(
            data = filter(data, is_selected),
            aes(label = gene_name),
            box.padding = 0.5,
            point.padding = 0.5,
            size = 4,
            color = "black",
            fontface = "bold.italic"
          )
        }
      } +
      scale_color_manual(
        values = c(
          "Upregulated" = "red",
          "Downregulated" = "blue",
          "Not Significant" = "gray"
        )
      ) +
      theme_minimal() +
      labs(
        title = paste(
          gsub("_", " ", input$volcano_comparison),
          "-",
          switch(
            input$volcano_sex,
            "both" = "Both sexes",
            "male" = "Male",
            "female" = "Female"
          )
        ),
        x = "log2 Fold Change",
        y = "-log10(P-Value)"
      ) +
      theme(legend.position = "none") +
      xlim(min_x - 0.5, max_x + 0.5)
    
    ggplotly(p, tooltip = "text") %>%
      layout(
        xaxis = list(title = "log2 Fold Change"),
        yaxis = list(title = "-log10(P-Value)")
      )
  })
  
  # -----------------------------------------------------
  # DEG tables
  # -----------------------------------------------------
  filtered_DEGs <- reactive({
    data <- data_objects()$DEGs_table
    
    if (input$comparison_filter != "All") {
      data <- data %>% filter(comparison == input$comparison_filter)
    }
    
    if (input$table_sex_filter != "All") {
      data <- data %>% filter(sex == input$table_sex_filter)
    }
    
    if (input$gene_name_filter != "All") {
      data <- data %>% filter(gene_name == input$gene_name_filter)
    }
    
    if (input$analysis_filter != "All") {
      data <- data %>% filter(analysis == input$analysis_filter)
    }
    
    keep_cols <- c("analysis", "sex", "comparison", "gene_id", "gene_name", "logFC", "AveExpr", "P.Value", "adj.P.Val")
    keep_cols <- intersect(keep_cols, colnames(data))
    
    data %>% select(all_of(keep_cols))
  })
  
  output$table <- DT::renderDataTable({
    DT::datatable(
      filtered_DEGs(),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      formatStyle("logFC", backgroundColor = styleInterval(brks, clrs)) %>%
      formatStyle("adj.P.Val", backgroundColor = styleInterval(brks2, clrs2)) %>%
      formatRound(columns = intersect(c("logFC", "AveExpr"), colnames(filtered_DEGs())), digits = 3) %>%
      formatRound(columns = intersect(c("P.Value", "adj.P.Val"), colnames(filtered_DEGs())), digits = 4)
  })
  
  # -----------------------------------------------------
  # CPM counts plots
  # -----------------------------------------------------
  output$counts_overall <- renderPlotly({
    req(input$counts_labels)
    
    gene_name <- input$counts_labels
    dat <- load_gene_counts(
      gene_name = gene_name,
      metadata = data_objects()$metadata,
      cpm_full = data_objects()$cpm_full_all,
      sex_subset = NULL
    )
    
    make_counts_subplot(dat, gene_name, title_prefix = "All samples")
  })
  
  output$counts_male <- renderPlotly({
    req(input$counts_labels)
    
    gene_name <- input$counts_labels
    dat <- load_gene_counts(
      gene_name = gene_name,
      metadata = data_objects()$metadata,
      cpm_full = data_objects()$cpm_full_male,
      sex_subset = "male"
    )
    
    make_counts_subplot(dat, gene_name, title_prefix = "Male")
  })
  
  output$counts_female <- renderPlotly({
    req(input$counts_labels)
    
    gene_name <- input$counts_labels
    dat <- load_gene_counts(
      gene_name = gene_name,
      metadata = data_objects()$metadata,
      cpm_full = data_objects()$cpm_full_female,
      sex_subset = "female"
    )
    
    make_counts_subplot(dat, gene_name, title_prefix = "Female")
  })
}

# ---------------------------------------------------------
# 5. RUN APP
# ---------------------------------------------------------
shinyApp(ui = ui, server = server)