# =========================================================================
# CROSS-PLATFORM SYSTEM INTEGRATION: COMPILING THE MASTER STUDY YIELD
# Layout: Unified Data Harmonization Layer
# =========================================================================
library(dplyr)
library(tidyr)
library(openxlsx) 
library(data.table)
library(pheatmap)

cat("--- Compiling Cross-Platform Global Yield Summaries ---\n")
base_dir <- "~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metabolic_Project/"

# Safely load analytical channel summary matrices
limma_summary_table_hilicneg <- fread(paste0(base_dir, "Hilic_Neg/Limma_summary_table.txt"))
limma_summary_table_hilicpos <- fread(paste0(base_dir, "Hilic_Pos/Limma_summary_table.txt"))
limma_summary_table_rpneg    <- fread(paste0(base_dir, "RP_Neg/Limma_summary_table.txt"))
limma_summary_table_rppos    <- fread(paste0(base_dir, "RP_Pos/Limma_summary_table.txt"))

# Add unified identifier tags to each platform data frame
yield_hilic_neg <- limma_summary_table_hilicneg %>% mutate(Platform = "HILIC_Negative")
yield_hilic_pos <- limma_summary_table_hilicpos %>% mutate(Platform = "HILIC_Positive") 
yield_rp_neg    <- limma_summary_table_rpneg %>% mutate(Platform = "RP_Negative")
yield_rp_pos    <- limma_summary_table_rppos %>% mutate(Platform = "RP_Positive")

# Combine datasets into a single master study structure
master_study_yield <- rbind(yield_hilic_neg, yield_hilic_pos, yield_rp_neg, yield_rp_pos)

# Pivot data into a wide, publication-ready summary grid
formatted_table <- master_study_yield %>%
  mutate(Summary_Text = paste0("FDR: ", Features_Strict_FDR, " | p<0.01: ", Features_Raw_P_01)) %>%
  select(Clinical_Trait, Platform, Summary_Text) %>%
  tidyr::pivot_wider(names_from = Platform, values_from = Summary_Text)

print("MASTER MANUSCRIPT MATRIX YIELD TABLE:")
print(formatted_table)

# Export the matrix directly to your project OneDrive directory
write.table(formatted_table, paste0(base_dir, "Master_Study_Yield_Matrix.txt"), sep="\t", row.names=FALSE, quote=FALSE)


# =========================================================================
# MULTI-PLATFORM INTEGRATED PHENOTYPIC CLUSTERING ENGINE
# Loop System: Full Cohort Clinical Trait Automation (A4 Landscape Vector)
# =========================================================================
cat("\n--- Initiating Multi-Trait Automation Engine ---\n")

# Core Sanitizer: Strips spaces, punctuation, and casing to ensure flawless lookups
sanitize   <- function(x) gsub("[[:space:][:punct:]_]", "", tolower(as.character(x)))
clean_text <- function(x) iconv(as.character(x), to = "UTF-8", sub = " ")

# Load Clinical Metadata
meta_stats <- fread(paste0(base_dir, "Metadata.txt"), sep = "\t", data.table = FALSE)
meta_stats$Sample_code <- trimws(as.character(meta_stats$Sample_code))
rownames(meta_stats)   <- meta_stats$Sample_code

# Standardized list of all analytical traits to loop through
target_traits <- c("Group", "Sex", "Disease_activity", "Disease_Type", "ANA", 
                   "HLA_B27", "Iritis", "Disease_outcome", "TMJ_arthritis", "Cervical_arthritis")

# Helper function: Scale features (rows) while keeping Samples in Columns
scale_platform_matrix <- function(mat_raw) {
  mtx_log <- log2(as.matrix(mat_raw) + 1)
  mtx_scaled <- t(apply(mtx_log, 1, function(x) {
    row_sd <- sd(x, na.rm = TRUE)
    if (is.na(row_sd) || row_sd == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / sqrt(row_sd)
  }))
  colnames(mtx_scaled) <- colnames(mat_raw)
  return(as.data.frame(mtx_scaled))
}

expression_matrices <- list()
feature_name_maps   <- list()

# Standardized Platform Files Map with underscore names to match Step 1 exactly
platform_files <- list(
  "HILIC Negative" = "Hilic_Neg/3_Clean_Data.txt",
  "HILIC Positive" = "Hilic_Pos/3_Clean_Data.txt",
  "RP Negative"    = "RP_Neg/3_Clean_Data.txt",
  "RP Positive"    = "RP_Pos/3_Clean_Data.txt"
)

# Ingest raw datasets and resolve sample code formatting variations dynamically
for (plat_name in names(platform_files)) {
  path <- paste0(base_dir, platform_files[[plat_name]])
  if (file.exists(path)) {
    df <- fread(path, sep = "\t", quote = "", header = TRUE, check.names = FALSE)
    
    clean_meta_codes <- sanitize(meta_stats$Sample_code)
    clean_df_headers <- sanitize(colnames(df))
    patient_indices  <- which(clean_df_headers %in% clean_meta_codes)
    
    if (length(patient_indices) > 0) {
      mtx <- as.matrix(df[, patient_indices, with = FALSE])
      colnames(mtx) <- meta_stats$Sample_code[match(clean_df_headers[patient_indices], clean_meta_codes)]
      
      final_id  <- clean_text(df[[1]])
      chem_name <- clean_text(df[[2]])
      chem_name[is.na(chem_name) | chem_name == ""] <- final_id[is.na(chem_name) | chem_name == ""]
      rownames(mtx) <- final_id
      
      expression_matrices[[plat_name]] <- scale_platform_matrix(mtx)
      feature_name_maps[[plat_name]]   <- data.frame(ID = final_id, Name = chem_name, stringsAsFactors = FALSE)
      cat(paste0("-> Loaded platform: ", plat_name, " (", ncol(mtx), " aligned samples)\n"))
    }
  }
}

common_samples <- Reduce(intersect, lapply(expression_matrices, colnames))
if (length(common_samples) < 2) stop("[CRITICAL CRASH]: Intersection profile holds insufficient sample depth.")

# Text-wrapping formatting tool for chemical names
wrap_chem_names <- function(names, width = 38) {
  sapply(names, function(x) {
    spaced <- gsub("-", "- ", x)
    wrapped <- paste(strwrap(spaced, width = width), collapse = "\n")
    cleaned <- gsub("- \n", "-\n", wrapped)
    cleaned <- gsub("- ", "-", cleaned)
    return(cleaned)
  })
}

expr_colors <- colorRampPalette(c("#0A2540", "#F4F6F8", "#B7094C"))(100)

# =========================================================================
# CORE EXECUTION LOOP: AUTOMATED TRAIT VISUALIZATION MATRIX
# =========================================================================
for (focus_trait in target_traits) {
  cat(paste0("\nProcessing Analytical Pipeline for Trait: [", focus_trait, "]\n"))
  
  master_inventory_list <- list()
  for (p in names(platform_files)) {
    stats_file <- paste0(base_dir, gsub("3_Clean_Data.txt", "", platform_files[[p]]), "Stats_All_Metabolites_", focus_trait, ".txt")
    if (file.exists(stats_file)) {
      stats <- fread(stats_file, sep = "\t", quote = "", fill = TRUE, check.names = FALSE)
      if (nrow(stats) > 0) {
        master_inventory_list[[p]] <- stats %>% filter(P.Value < 0.05) %>% mutate(Platform = p)
      }
    }
  }
  
  if (length(master_inventory_list) == 0) {
    cat(paste0("-> Skipping [", focus_trait, "]: No statistical discovery tables found.\n"))
    next
  }
  
  master_inventory <- rbindlist(master_inventory_list, fill = TRUE) %>% arrange(P.Value) %>% head(40)
  if (nrow(master_inventory) == 0) {
    cat(paste0("-> Skipping [", focus_trait, "]: Statistical inventory holds 0 significant metrics.\n"))
    next
  }
  
  heatmap_expr_list <- list()
  feature_platforms <- c()
  
  for (i in 1:nrow(master_inventory)) {
    f_id   <- as.character(master_inventory$Metabolite_ID[i])
    p_name <- master_inventory$Platform[i]
    df     <- expression_matrices[[p_name]]
    catalog <- feature_name_maps[[p_name]]
    
    if (is.null(df) || is.null(catalog)) next
    
    row_idx <- which(catalog$ID == f_id | catalog$Name == f_id)[1]
    
    if (!is.na(row_idx)) {
      actual_id <- catalog$ID[row_idx]
      display_name <- catalog$Name[row_idx]
      unique_feature_name <- paste0("[", p_name, "] ", display_name)
      
      heatmap_expr_list[[unique_feature_name]] <- as.numeric(df[actual_id, common_samples])
      feature_platforms <- c(feature_platforms, p_name)
    }
  }
  
  if (length(heatmap_expr_list) == 0) {
    cat(paste0("-> Warning: Feature cross-referencing returned empty matrices for trait: ", focus_trait, "\n"))
    next
  }
  
  heatmap_matrix <- do.call(rbind, heatmap_expr_list)
  colnames(heatmap_matrix) <- common_samples
  
  trait_vector <- meta_stats[common_samples, focus_trait]
  valid_canvas_patients <- common_samples[!is.na(trait_vector) & trait_vector != ""]
  
  if (length(valid_canvas_patients) < 2) {
    cat(paste0("-> Skipping [", focus_trait, "]: Insufficient active patients without NA fields.\n"))
    next
  }
  
  heatmap_matrix_sub <- heatmap_matrix[, valid_canvas_patients, drop = FALSE]
  meta_stats_sub     <- meta_stats[valid_canvas_patients, ]
  
  # Structural Safety Lock: Impute missing entries with the relative cohort mean (0)
  heatmap_matrix_sub[is.na(heatmap_matrix_sub) | is.nan(heatmap_matrix_sub) | is.infinite(heatmap_matrix_sub)] <- 0
  
  # Adaptive Variance Filter: Exclude invariant features inside small sub-cohort cuts
  row_variances <- apply(heatmap_matrix_sub, 1, sd, na.rm = TRUE)
  active_rows   <- which(!is.na(row_variances) & row_variances > 1e-4)
  
  if (length(active_rows) < 2) {
    cat(paste0("-> Skipping [", focus_trait, "]: Less than 2 active variations inside sub-cohort.\n"))
    next
  }
  
  heatmap_matrix_sub <- heatmap_matrix_sub[active_rows, , drop = FALSE]
  feature_platforms  <- feature_platforms[active_rows]
  
  # Standardize primary factor levels securely using direct construction injection
  if (focus_trait == "Disease_Type") {
    primary_factor <- factor(meta_stats_sub[[focus_trait]], levels = c(0, 1, 2, 3), labels = c("SLA", "ERA", "Oligo", "Poly"))
  } else if (focus_trait == "Disease_activity") {
    primary_factor <- factor(meta_stats_sub[[focus_trait]], levels = c(1, 2, 3), labels = c("Low", "Moderate", "High"))
  } else if (focus_trait == "Sex") {
    primary_factor <- factor(meta_stats_sub[[focus_trait]], levels = c(1, 2), labels = c("Female", "Male"))
  } else {
    primary_factor <- factor(meta_stats_sub[[focus_trait]])
  }
  
  # Secure data frame initialization with the first column bound directly
  ann_col <- data.frame(Target_Trait = primary_factor, row.names = valid_canvas_patients)
  colnames(ann_col) <- focus_trait
  
  # Append context tracking tracks cleanly if they aren't the primary variable
  if (focus_trait != "Sex" && "Sex" %in% colnames(meta_stats_sub)) {
    ann_col$Sex <- factor(meta_stats_sub$Sex, levels = c(1, 2), labels = c("Female", "Male"))
  }
  if (focus_trait != "Disease_activity" && "Disease_activity" %in% colnames(meta_stats_sub)) {
    ann_col$Activity <- factor(meta_stats_sub$Disease_activity, levels = c(1, 2, 3), labels = c("Low", "Moderate", "High"))
  }
  
  ann_row <- data.frame(Platform = factor(feature_platforms), row.names = rownames(heatmap_matrix_sub))
  
  # Execute manual pre-sorting hierarchical clustering coordinates
  dist_cols <- dist(t(heatmap_matrix_sub), method = "euclidean")
  hc_cols   <- hclust(dist_cols, method = "ward.D2")
  cluster_order <- hc_cols$order
  
  heatmap_final <- heatmap_matrix_sub[, cluster_order, drop = FALSE]
  ann_final     <- ann_col[cluster_order, , drop = FALSE]
  row_display_labels <- wrap_chem_names(rownames(heatmap_final))
  
  save_filename <- paste0(base_dir, "Automated_Heatmap_Trait_", focus_trait, "_A4.pdf")
  
  pheatmap(
    mat = as.matrix(heatmap_final),
    annotation_col = ann_final,
    annotation_row = ann_row,
    cluster_cols = FALSE,               # Order is pre-sorted manually
    cluster_rows = TRUE,                
    color = expr_colors,
    border_color = NA,
    labels_row = row_display_labels,
    fontsize_row = 5.2,
    filename = save_filename,
    width = 11.69, height = 8.27        
  )
  cat(paste0("SUCCESS: Exported vector document at: ", save_filename, "\n"))
}

cat("\n=========================================================================\n")
cat("TRAIT AUTOMATION COMPLETE. PORTFOLIO OF DOCUMENTATION COMPILED SUCCESSFULLY.\n")
cat("=========================================================================\n")
