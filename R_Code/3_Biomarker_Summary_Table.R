# =========================================================================
# MANUSCRIPT ARCHITECTURE: COMPILING THE MASTER BIOMARKER SUMMARY TABLE
# Layout: Standard Publication-Grade Results Matrix (Tab-Delimited Export)
# =========================================================================
library(dplyr)
library(tidyr)
library(data.table)

cat("--- Initiating Master Biomarker Extraction Engine ---\n")
base_dir <- "~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metabolic_Project/"

# Core Sanitizer: Ensures clean character translations for chemical entities
clean_text <- function(x) iconv(as.character(x), to = "UTF-8", sub = " ")
sanitize   <- function(x) gsub("[[:space:][:punct:]_]", "", tolower(as.character(x)))

# Standardized list of traits that successfully yielded data in your cohort
target_traits <- c("Sex", "Disease_activity", "Disease_Type", "ANA", "HLA_B27", "TMJ_arthritis")

# -------------------------------------------------------------------------
# Stage 1: Build Global Feature Chemical Maps
# -------------------------------------------------------------------------
cat("Building cross-platform chemical nomenclature dictionaries...\n")
feature_name_maps <- list()
platform_files <- list(
  "HILIC_Negative" = "Hilic_Neg/3_Clean_Data.txt",
  "HILIC_Positive" = "Hilic_Pos/3_Clean_Data.txt",
  "RP_Negative"    = "RP_Neg/3_Clean_Data.txt",
  "RP_Positive"    = "RP_Pos/3_Clean_Data.txt"
)

for (plat_name in names(platform_files)) {
  path <- paste0(base_dir, platform_files[[plat_name]])
  if (file.exists(path)) {
    df <- fread(path, sep = "\t", quote = "", header = TRUE, select = c(1,2), check.names = FALSE)
    final_id  <- clean_text(df[[1]])
    chem_name <- clean_text(df[[2]])
    # Fallback to ID if chemical name column is missing or blank
    chem_name[is.na(chem_name) | chem_name == ""] <- final_id[is.na(chem_name) | chem_name == ""]
    
    feature_name_maps[[plat_name]] <- data.frame(ID = final_id, Name = chem_name, stringsAsFactors = FALSE)
  }
}

# -------------------------------------------------------------------------
# Stage 2: Iterate and Extract Top Significant Entities per Trait
# -------------------------------------------------------------------------
master_table_list <- list()

for (focus_trait in target_traits) {
  cat(paste0("Extracting statistical vectors for Trait: [", focus_trait, "]\n"))
  
  trait_hits_list <- list()
  for (p in names(platform_files)) {
    stats_file <- paste0(base_dir, gsub("3_Clean_Data.txt", "", platform_files[[p]]), "Stats_All_Metabolites_", focus_trait, ".txt")
    
    if (file.exists(stats_file)) {
      stats <- fread(stats_file, sep = "\t", quote = "", fill = TRUE, check.names = FALSE)
      if (nrow(stats) > 0) {
        # Enforce strict significance filtering
        sig_stats <- stats %>% filter(P.Value < 0.05) %>% mutate(Platform = p, Clinical_Trait = focus_trait)
        trait_hits_list[[p]] <- sig_stats
      }
    }
  }
  
  if (length(trait_hits_list) == 0) next
  
  # Bind all platforms for this specific trait and sort by raw significance
  trait_master <- rbindlist(trait_hits_list, fill = TRUE) %>% arrange(P.Value)
  
  # Cap at the top 10 standout biomarkers per trait to keep the paper table concise
  trait_top10 <- head(trait_master, 10)
  
  if (nrow(trait_top10) == 0) next
  
  # Resolve IDs to clean chemical designations
  resolved_rows <- list()
  for (i in 1:nrow(trait_top10)) {
    f_id    <- as.character(trait_top10$Metabolite_ID[i])
    p_name  <- trait_top10$Platform[i]
    catalog <- feature_name_maps[[p_name]]
    
    # Locate row index via Dual-Key Lookup mapping
    row_idx <- which(catalog$ID == f_id | catalog$Name == f_id)[1]
    display_name <- if(!is.na(row_idx)) catalog$Name[row_idx] else f_id
    
    # Dynamically extract effect size / direction columns depending on limma configurations
    logFC_val <- if("logFC" %in% colnames(trait_top10)) trait_top10$logFC[i] else NA
    t_stat    <- if("t" %in% colnames(trait_top10)) trait_top10$t[i] else NA
    adj_P     <- if("adj.P.Val" %in% colnames(trait_top10)) trait_top10$adj.P.Val[i] else NA
    
    resolved_rows[[i]] <- data.frame(
      Clinical_Trait    = focus_trait,
      Metabolite_Name   = display_name,
      Platform          = p_name,
      Log_Fold_Change   = logFC_val,
      T_Statistic       = t_stat,
      Raw_P_Value       = trait_top10$P.Value[i],
      Adjusted_FDR      = adj_P,
      stringsAsFactors  = FALSE
    )
  }
  
  master_table_list[[focus_trait]] <- rbindlist(resolved_rows)
}

# -------------------------------------------------------------------------
# Stage 3: Final Consolidation and Aesthetic Formatting
# -------------------------------------------------------------------------
final_summary_table <- rbindlist(master_table_list)

# Format numeric strings to meet rigorous style guidelines for elite journals
final_summary_table <- final_summary_table %>%
  mutate(
    # 1. Force strict numeric coercion into intermediate tracking variables
    # This strips out any hidden character types causing the sprintf crash
    Num_LogFC = as.numeric(as.character(Log_Fold_Change)),
    Num_TStat = as.numeric(as.character(T_Statistic)),
    Num_PVal  = as.numeric(as.character(Raw_P_Value)),
    Num_FDR   = as.numeric(as.character(Adjusted_FDR))
  ) %>%
  mutate(
    # 2. Safely apply scientific formatting on the guaranteed numeric variables
    Log_Fold_Change = ifelse(is.na(Num_LogFC), "N/A", sprintf("%.3f", Num_LogFC)),
    T_Statistic     = ifelse(is.na(Num_TStat), "N/A", sprintf("%.2f", Num_TStat)),
    Raw_P_Value     = ifelse(is.na(Num_PVal),  "N/A", sprintf("%.2e", Num_PVal)),
    Adjusted_FDR    = ifelse(is.na(Num_FDR),   "N/A", sprintf("%.2e", Num_FDR))
  ) %>%
  # 3. Clean up the workspace by dropping the intermediate helper columns
  select(-Num_LogFC, -Num_TStat, -Num_PVal, -Num_FDR)

# Export tab-separated spreadsheet file
output_path <- paste0(base_dir, "Master_Biomarker_Summary_Table.txt")
write.table(final_summary_table, output_path, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=========================================================================\n")
cat("MASTER BIOMARKER TABLE COMPILED SUCCESSFULLY!\n")
cat(paste0("Location: ", output_path, "\n"))
cat("=========================================================================\n")
