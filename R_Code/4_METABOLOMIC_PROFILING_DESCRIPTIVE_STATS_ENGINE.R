# =========================================================================
# MANUSCRIPT ARCHITECTURE: METABOLOMIC PROFILING & DESCRIPTIVE STATS ENGINE
# Layout: Unified Multi-Trait Boxplot Panels & Synchronized Data Tables
# =========================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(gridExtra)

cat("--- Initiating Optimized Global Ingestion & Visual Profiling Engine ---\n")
base_dir <- "~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metabolic_Project/"

# Core Sanitizer tools
clean_text <- function(x) iconv(as.character(x), to = "UTF-8", sub = " ")
sanitize   <- function(x) gsub("[[:space:][:punct:]_]", "", tolower(as.character(x)))

# Load Metadata using fault-tolerant stream parser
meta_stats <- read.delim(paste0(base_dir, "Metadata.txt"), sep = "\t", header = TRUE, check.names = FALSE, fill = TRUE, quote = "")
meta_stats$Sample_code <- trimws(as.character(meta_stats$Sample_code))
rownames(meta_stats)   <- meta_stats$Sample_code

# -------------------------------------------------------------------------
# Stage 1: Ingest and Normalize Abundances (Fault-Tolerant Matrix Lock)
# -------------------------------------------------------------------------
cat("Loading cross-platform data channels into memory...\n")
expression_matrices <- list()
feature_name_maps   <- list()

platform_files <- list(
  "HILIC_Negative" = "Hilic_Neg/3_Clean_Data.txt",
  "HILIC_Positive" = "Hilic_Pos/3_Clean_Data.txt",
  "RP_Negative"    = "RP_Neg/3_Clean_Data.txt",
  "RP_Positive"    = "RP_Pos/3_Clean_Data.txt"
)

for (plat_name in names(platform_files)) {
  path <- paste0(base_dir, platform_files[[plat_name]])
  if (file.exists(path)) {
    df <- read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, fill = TRUE, quote = "")
    
    clean_meta_codes <- sanitize(meta_stats$Sample_code)
    clean_df_headers <- sanitize(colnames(df))
    patient_indices  <- which(clean_df_headers %in% clean_meta_codes)
    
    if (length(patient_indices) > 0) {
      mtx <- as.matrix(df[, patient_indices, drop = FALSE])
      colnames(mtx) <- meta_stats$Sample_code[match(clean_df_headers[patient_indices], clean_meta_codes)]
      
      final_id  <- clean_text(df[[1]])
      chem_name <- clean_text(df[[2]])
      chem_name[is.na(chem_name) | chem_name == ""] <- final_id[is.na(chem_name) | chem_name == ""]
      
      valid_rows <- !is.na(final_id) & final_id != ""
      mtx        <- mtx[valid_rows, , drop = FALSE]
      final_id   <- final_id[valid_rows]
      chem_name  <- chem_name[valid_rows]
      
      numeric_mtx <- matrix(as.numeric(mtx), nrow = nrow(mtx), ncol = ncol(mtx))
      colnames(numeric_mtx) <- colnames(mtx)
      rownames(numeric_mtx) <- final_id
      
      row_has_data <- rowSums(!is.na(numeric_mtx)) > 0
      numeric_mtx  <- numeric_mtx[row_has_data, , drop = FALSE]
      final_id     <- final_id[row_has_data]
      chem_name    <- chem_name[row_has_data]
      
      expression_matrices[[plat_name]] <- log2(numeric_mtx + 1)
      feature_name_maps[[plat_name]]   <- data.frame(ID = final_id, Name = chem_name, stringsAsFactors = FALSE)
      cat(paste0("   -> Platform [", plat_name, "] cached: ", nrow(numeric_mtx), " rows passed.\n"))
    }
  }
}

common_samples <- Reduce(intersect, lapply(expression_matrices, colnames))
cat(paste0("STATUS: Found [", length(common_samples), "] fully aligned patient profiles.\n\n"))

# -------------------------------------------------------------------------
# Define Optimized High-Contrast Chromatic Palettes for Presentation
# -------------------------------------------------------------------------
color_palettes <- list(
  "Sex"              = c("#E11D48", "#0891B2"),            # Premium Deep Coral / Deep Cyan
  "Disease_activity" = c("#10B981", "#0284C7", "#DC2626"), # Emerald / Deep Oceanic Blue / Vivid Crimson
  "Disease_Type"     = c("#475569", "#1E1B4B", "#6D28D9", "#C084FC"), # Slate / Indigo / Royal Violet / Orchid
  "ANA"              = c("#4B5563", "#D97706"),            # Charcoal Muted Grey / Polished Amber Gold
  "HLA_B27"          = c("#64748B", "#2563EB"),            # Cool Steel Grey / High-Contrast Electric Blue
  "TMJ_arthritis"    = c("#65A30D", "#BE123C")             # Clear Olive Green / Dark Ruby Crimson
)

# Global accumulator for descriptive stats (Filled during unified execution loop)
master_descriptive_dump <- list()
all_traits <- c("Sex", "Disease_activity", "Disease_Type", "ANA", "HLA_B27", "TMJ_arthritis")

for (focus_trait in all_traits) {
  cat(paste0("=================== PROCESSING TRAIT: ", focus_trait, " ===================\n"))
  
  # Stage 2: Identify Significant Features
  master_inventory_list <- list()
  for (p in names(platform_files)) {
    stats_file <- paste0(base_dir, gsub("3_Clean_Data.txt", "", platform_files[[p]]), "Stats_All_Metabolites_", focus_trait, ".txt")
    if (file.exists(stats_file)) {
      stats <- read.delim(stats_file, sep = "\t", header = TRUE, check.names = FALSE, fill = TRUE, quote = "")
      if (nrow(stats) > 0 && "P.Value" %in% colnames(stats)) {
        stats$P.Value <- as.numeric(as.character(stats$P.Value))
        master_inventory_list[[p]] <- stats %>% filter(P.Value < 0.05) %>% mutate(Platform = p)
      }
    }
  }
  
  if (length(master_inventory_list) == 0) {
    cat(paste0("   ! Skipping trait [", focus_trait, "]: No statistical input text tables found.\n"))
    next
  }
  
  master_inventory <- rbindlist(master_inventory_list, fill = TRUE) %>% arrange(P.Value) %>% head(6)
  if (nrow(master_inventory) == 0) {
    cat(paste0("   ! Skipping trait [", focus_trait, "]: Zero molecules passed the alpha threshold.\n"))
    next
  }
  
  # Stage 3: Extract Data Points & Apply Chemical Title Alignment Engine
  plot_data_list <- list()
  
  for (i in 1:nrow(master_inventory)) {
    f_id    <- as.character(master_inventory$Metabolite_ID[i])
    p_name  <- master_inventory$Platform[i]
    p_val   <- master_inventory$P.Value[i]
    df      <- expression_matrices[[p_name]] 
    catalog <- feature_name_maps[[p_name]]
    
    if (is.null(df) || is.null(catalog)) next
    
    row_idx <- which(catalog$ID == f_id | catalog$Name == f_id | tolower(catalog$ID) == tolower(f_id) | tolower(catalog$Name) == tolower(f_id))[1]
    if (is.na(row_idx)) {
      clean_f_id <- sanitize(f_id)
      row_idx <- which(sanitize(catalog$ID) == clean_f_id | sanitize(catalog$Name) == clean_f_id)[1]
    }
    if (is.na(row_idx)) {
      row_idx <- grep(f_id, catalog$Name, fixed = TRUE)[1]
      if (is.na(row_idx)) row_idx <- grep(f_id, catalog$ID, fixed = TRUE)[1]
    }
    
    if (!is.na(row_idx)) {
      actual_id    <- catalog$ID[row_idx]
      display_name <- catalog$Name[row_idx]
      
      # SMART CHEMICAL TITLE ENGINE (Sanitizes names for both Boxplots AND Text Summaries)
      display_name <- gsub('^"|"$', '', display_name)
      display_name <- gsub("^'|'$", "", display_name)
      display_name <- trimws(display_name)
      
      if (startsWith(display_name, ")") || startsWith(display_name, "]")) {
        display_name <- substring(display_name, 2)
      }
      if (nchar(display_name) > 42) {
        display_name <- paste0(substr(display_name, 1, 39), "...")
      }
      
      wrapped_name <- paste(strwrap(display_name, width = 32), collapse = "\n")
      plot_title   <- paste0("[", gsub("_", " ", p_name), "]\n", wrapped_name, "\np = ", sprintf("%.2e", p_val))
      
      abundance <- df[actual_id, common_samples]
      trait_raw <- meta_stats[common_samples, focus_trait]
      
      feature_df <- data.frame(
        Sample       = common_samples,
        Abundance    = as.numeric(abundance),
        Trait_Group  = trimws(as.character(trait_raw)), 
        Feature_Key  = plot_title,
        stringsAsFactors = FALSE
      ) %>% filter(!is.na(Trait_Group) & Trait_Group != "" & Trait_Group != "NA")
      
      if (nrow(feature_df) < 2) next
      
      # Harmonize levels and assign clear categorical metadata definitions
      if (focus_trait == "Disease_activity") {
        feature_df$Trait_Group <- factor(feature_df$Trait_Group, levels = c("1", "2", "3"), labels = c("Low", "Moderate", "High"))
      } else if (focus_trait == "Disease_Type") {
        feature_df$Trait_Group <- factor(feature_df$Trait_Group, levels = c("0", "1", "2", "3"), labels = c("SLA", "ERA", "Oligo", "Poly"))
      } else if (focus_trait == "Sex") {
        feature_df$Trait_Group <- factor(feature_df$Trait_Group, levels = c("1", "2"), labels = c("Female", "Male"))
      } else if (focus_trait %in% c("ANA", "HLA_B27", "TMJ_arthritis")) {
        feature_df$Trait_Group <- factor(feature_df$Trait_Group, levels = c("0", "1"), labels = c("Negative", "Positive"))
      } else {
        feature_df$Trait_Group <- factor(feature_df$Trait_Group)
      }
      
      # CRITICAL ACCELERATION LOCK: Compute descriptive metrics here instantly
      summary_metrics <- feature_df %>%
        group_by(Trait_Group) %>%
        summarise(
          Cohort_N     = n(),
          Mean_Log2    = round(mean(Abundance, na.rm = TRUE), 4),
          SD_Log2      = round(sd(Abundance, na.rm = TRUE), 4),
          Median_Log2  = round(median(Abundance, na.rm = TRUE), 4),
          IQR_Log2     = round(IQR(Abundance, na.rm = TRUE), 4),
          Min_Log2     = round(min(Abundance, na.rm = TRUE), 4),
          Max_Log2     = round(max(Abundance, na.rm = TRUE), 4),
          .groups      = 'drop'
        ) %>%
        mutate(
          Clinical_Trait  = focus_trait,
          Metabolite_ID   = actual_id,
          Chemical_Name   = display_name,   # 100% matched to your boxplot text limits
          Platform        = gsub("_", " ", p_name),
          ANOVA_P_Value   = sprintf("%.2e", p_val)
        ) %>%
        select(Clinical_Trait, Metabolite_ID, Chemical_Name, Platform, ANOVA_P_Value, everything())
      
      master_descriptive_dump[[length(master_descriptive_dump) + 1]] <- summary_metrics
      plot_data_list[[length(plot_data_list) + 1]] <- feature_df
    }
  }
  
  # Stage 4: Render High-Contrast Graphics (With Correct Else Block Lock)
  if (length(plot_data_list) == 0) {
    cat("   ! Target feature profiles returned empty sets.\n")
  } else {
    boxplot_objects <- list()
    active_palette  <- color_palettes[[focus_trait]]
    
    for (j in seq_along(plot_data_list)) {
      df_sub       <- plot_data_list[[j]]
      title_string <- df_sub$Feature_Key[1]
      
      p <- ggplot(df_sub, aes(x = Trait_Group, y = Abundance, fill = Trait_Group)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.45, color = "#1E293B", lwd = 0.5) +
        geom_jitter(width = 0.12, size = 1.4, alpha = 0.55, aes(color = Trait_Group)) +
        labs(title = title_string, x = NULL, y = "Relative Abundance (Log2)") +
        theme_classic(base_size = 9) +
        theme(
          plot.title         = element_text(face = "bold", size = 7.5, hjust = 0.5, lineheight = 1.05, color = "#0F172A"),
          axis.text.x        = element_text(face = "bold", color = "#334155", size = 8),
          axis.text.y        = element_text(color = "#475569"),
          axis.title.y       = element_text(size = 8, color = "#475569"),
          legend.position    = "none",
          panel.grid.major.y = element_line(color = "#F1F5F9", linewidth = 0.4),
          plot.margin        = margin(10, 5, 10, 5)
        )
      
      if (!is.null(active_palette)) {
        p <- p + scale_fill_manual(values = active_palette) + scale_color_manual(values = active_palette)
      }
      boxplot_objects[[j]] <- p
    }
    
    total_plots    <- length(boxplot_objects)
    columns_needed <- min(total_plots, 3)
    rows_needed    <- ceiling(total_plots / columns_needed)
    
    output_path <- paste0(base_dir, "Top_Biomarker_Profiles_", focus_trait, ".pdf")
    pdf(output_path, width = 11.69, height = 8.27) 
    grid.arrange(grobs = boxplot_objects, ncol = columns_needed, nrow = rows_needed)
    dev.off()
    
    cat(paste0("   -> SUCCESS: Compiled vector PDF panels -> ", output_path, "\n\n"))
  } # <-- FIX: This brace safely closes the Stage 4 'else' statement
} # <-- FIX: This brace safely terminates the trait processing block loop

# -------------------------------------------------------------------------
# Compile, Export and Finalize Structured Tables
# -------------------------------------------------------------------------
cat("--- Finalizing Supplementary Matrix Deliverables ---\n")
if (length(master_descriptive_dump) > 0) {
  final_descriptive_table <- rbindlist(master_descriptive_dump, fill = TRUE)
  output_table_path <- paste0(base_dir, "Top_Biomarkers_Descriptive_Group_Statistics.txt")
  write.table(final_descriptive_table, file = output_table_path, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("\n=========================================================================\n")
  cat("DESCRIPTIVE TABLES GENERATION COMPLETE!\n")
  cat(paste0("Summary Table Exported to: ", output_table_path, "\n"))
  cat("=========================================================================\n")
} else {
  cat("ERROR: Statistics cache empty. Check input data bounds.\n")
}
