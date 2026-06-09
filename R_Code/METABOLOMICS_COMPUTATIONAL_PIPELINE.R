# =========================================================================
# METABOLOMICS COMPUTATIONAL PIPELINE: RP Negative PLATFORM
# Workflow: Preprocessing, Global Screening (PERMANOVA) & Biomarker Discovery (Limma)
# =========================================================================

# -------------------------------------------------------------------------
# MODULE 1: ENVIRONMENT SETUP & DEPENDENCIES
# -------------------------------------------------------------------------
library(data.table)   # High-performance data reading and manipulation
library(vegan)        # Multivariate ecological/omics statistics (PERMANOVA)
library(ggplot2)      # Core publication-grade data visualization engine
library(viridis)      # Perceptually uniform, colorblind-friendly palettes
library(limma)        # Empirical Bayes moderated linear modeling for biomarker yield
library(dplyr)        # Elegant data frame plumbing and transformations
library(ggrepel)      # Smart text label placement avoiding graphic overlaps
library(tidyr)
library(patchwork) # Core package for combining multiple plots cleanly
# Suppress scientific notation globally for clean axis text formatting
options(scipen = 999)

# -------------------------------------------------------------------------
# MODULE 2: DATA ACQUISITION & PREPROCESSING PIPELINE
# -------------------------------------------------------------------------
# Load raw intensity feature matrix. 
# CRITICAL: check.names=FALSE protects raw chemical names with spaces/hyphens from scrambling.
df <- fread("~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metabolic_Project/RP_Neg/3_Clean_Data.txt", check.names = FALSE)

# Extract feature measurement zone (Samples are positioned across columns 4 to 65)
mtx <- df[, 4:68]

# Deduplicate original metabolite tags safely using make.unique
unique_names <- make.unique(as.character(df[[1]]))
rownames(mtx) <- unique_names

# --- Step A: Sample Normalization (The "Vertical" Correction) ---
# Aligns individual injection volumes to the global median to eliminate technical variation
sample_medians <- apply(mtx, 2, median, na.rm = TRUE)
mtx_norm <- sweep(mtx, 2, sample_medians, "/") * median(sample_medians)

# --- Step B: Data Transformation (The "Distribution" Correction) ---
# Log2 transformation dampens heteroscedasticity and converts highly skewed exponential ranges to normal distributions
# Adding a pseudo-count (+1) safely avoids mathematically undefined log(0) errors
mtx_log <- log2(mtx_norm + 1)

# --- Step C: Data Scaling (The "Horizontal" Correction) ---
# Pareto Scaling mean-centers rows and divides by the square root of the standard deviation.
# This prevents dominant high-abundance metabolites from drowning out low-abundance regulatory features.
mtx_final <- t(apply(mtx_log, 1, function(x) (x - mean(x)) / sqrt(sd(x))))

# --- Step D: Metadata Ingestion ---
metadata <- fread("~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metabolic_Project/Metadata_RP_Neg.txt")

# -------------------------------------------------------------------------
# MODULE 3: DIAGNOSTIC EXPLORATORY ANALYSIS (GLOBAL PCA)
# -------------------------------------------------------------------------
# Run Principal Component Analysis across unconstrained features
pca_res <- prcomp(t(mtx_final), center = FALSE, scale. = FALSE)

# Compute percentage of global variance captured by individual orthogonal axes
var_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

# Construct data coordinate frame mapping sample projections to clinical metadata
pca_df <- data.frame(SampleID = rownames(pca_res$x),
                     PC1 = pca_res$x[,1],
                     PC2 = pca_res$x[,2],
                     Group = metadata$Group)

# Render Diagnostic Score Plot
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA: Global Unconstrained Metabolic Profiling (RP Negative)",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  scale_color_brewer(palette = "Set1") + 
  stat_ellipse(level = 0.95) + # Draws 95% confidence intervals per cluster
  scale_x_continuous(breaks = round(seq(-50, 40, by = 5), 1)) +
  scale_y_continuous(breaks = round(seq(-25, 20, by = 5), 1)) +
  theme(plot.title = element_text(family = "serif", size = 18, face = "bold"),
        axis.title.x = element_text(family = "serif", size = 16),
        axis.title.y = element_text(family = "serif", size = 16),
        axis.text.x = element_text(family = "serif", size = 14),
        axis.text.y = element_text(family = "serif", size = 14),
        legend.title = element_text(family = "serif", size = 16),
        legend.text = element_text(family = "serif", size = 16),
        panel.background = element_blank())

print(p_pca)

# -------------------------------------------------------------------------
# MODULE 4: COHORT PURIFICATION & DATA ALIGNMENT
# -------------------------------------------------------------------------
# Transpose matrices so samples match row-wise. 
# CRITICAL: check.names=FALSE is passed again to protect special characters in column names.
norm_data <- as.data.frame(t(mtx_final), check.names = FALSE)

# Technical pool samples (QCs) are removed before clinical profiling to prevent variance inflation
# Ensure sequence in the meta file and data file is same
is_qc <- metadata$Group == "QC" 
mtx_stats <- norm_data[!is_qc, ]
meta_stats <- metadata[!is_qc, ]

# Register the 10 distinct clinical groupings to evaluate
groupings <- c("Age", "Sex", "Disease_activity", "Disease_Type", 
               "ANA", "HLA_B27", "Iritis", "Disease_outcome", 
               "TMJ_arthritis", "Cervical_arthritis")

# Programmatically cast clinical classifications to factors
for(g in groupings) {
  meta_stats[[g]] <- as.factor(meta_stats[[g]])
}

# -------------------------------------------------------------------------
# MODULE 5: PHASE I SCREENING — GLOBAL VARIANCE QUANTIFICATION (PERMANOVA)
# -------------------------------------------------------------------------
# Initialize container data frame to hold global separation scores
permanova_table <- data.frame(
  Clinical_Trait = character(),
  Variance_Explained_R2 = numeric(),
  P_Value = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

cat("--- Running Global PERMANOVA Screen ---\n")
for (g in groupings) {
  # Dynamically synchronize complete-cases across target factors, age, and sex covariates
  valid_samples <- !is.na(meta_stats[[g]]) & !is.na(meta_stats$Age) & !is.na(meta_stats$Sex)
  
  meta_sub <- meta_stats[valid_samples, ]
  mtx_sub <- mtx_stats[valid_samples, ]
  # Clean up factor levels from dropped groups
  meta_sub <- droplevels(meta_sub)
  
  # Safety skip to protect against clinical layers dropping below 2 comparative dimensions
  if (length(unique(na.omit(meta_stats[[g]]))) < 2) next
  
  # Execute non-parametric multivariate analysis of variance via Euclidean space distance
  set.seed(123) # Locks permuted data scrambling structure for exact reproducibility
  # STATISTICAL FIX: Multi-variable formula to control for Age and Sex globally
  if (g == "Age") 
    { formula_str <- "mtx_sub ~ Sex + Age" }
  else if (g == "Sex") 
    { formula_str <- "mtx_sub ~ Age + Sex" }
  else { formula_str <- paste("mtx_sub ~ Age + Sex +", g) }
  
  # Explicitly coerce the text string into a formal R formula object
  formula_object <- as.formula(formula_str)
  
  # by = "margin" isolates the independent effect size of the clinical trait
  res <- adonis2(formula_object, data = meta_sub, method = "euclidean", permutations = 999, by = "margin")
  
  # Extract true effect sizes (R2) and structural significance parameters
  trait_row <- grep(g, rownames(res))[1]
  r2 <- res$R2[trait_row]
  p_val <- res$`Pr(>F)`[trait_row]
  
  sig_label <- ifelse(p_val < 0.05, "* Significant", ifelse(p_val < 0.1, ". Trend", "Not Significant"))
  
  permanova_table <- rbind(permanova_table, data.frame(
    Clinical_Trait = g, Variance_Explained_R2 = r2, P_Value = p_val, Significance = sig_label
  ))
}

# Rank clinical parameters based on maximum total macro variance captured
permanova_table <- permanova_table[order(-permanova_table$Variance_Explained_R2), ]

print("PERMANOVA SUMMARY TABLE:")
print(permanova_table)
write.table (permanova_table, "~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metabolic_Project/RP_Neg/Permanova_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# --- Visualizing Phase I Summary (Continuous Viridis Ranking Plot) ---
# Arrange y-axis labels based on performance ranking hierarchy
permanova_table$Clinical_Trait <- factor(permanova_table$Clinical_Trait, levels = rev(permanova_table$Clinical_Trait))

p_perm <- ggplot(permanova_table, aes(x = Clinical_Trait, y = Variance_Explained_R2 * 100)) +
  geom_bar(aes(fill = Variance_Explained_R2 * 100), stat = "identity", color = "black", width = 0.65, alpha = 0.9) +
  geom_text(aes(label = paste0("p = ", sprintf("%.3f", P_Value))), 
            hjust = -0.2, family = "serif", size = 4.5, fontface = "italic", color = "grey20") +
  scale_fill_viridis_c(option = "plasma", name = "Global Variance\nExplained (%)") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Multivariate Screening: Ranking Patient Trait Influence (RP Negative)",
       subtitle = "Ranked by Global Variance Explained (PERMANOVA Euclidean Distance)",
       x = "Clinical Phenotype / Trait",
       y = "Total Percentage of Global Metabolic Variance Explained (%)") +
  theme(text = element_text(family = "serif"),
        plot.title = element_text(size = 16, face = "bold", color = "black"),
        plot.subtitle = element_text(size = 11, face = "italic", color = "grey40"),
        axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 12, face = "bold", color = "black"),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3)))

print(p_perm)

# -------------------------------------------------------------------------
# MODULE 6: PHASE II — COVARIATE-ADJUSTED BIOMARKER DISCOVERY (LIMMA)
# -------------------------------------------------------------------------
# Storage configurations for statistical arrays and numerical yields
all_limma_results <- list()
limma_summary_table <- data.frame(
  Clinical_Trait = character(),
  Features_Raw_P_05 = integer(),
  Features_Raw_P_01 = integer(),
  Features_Strict_FDR = integer(),
  stringsAsFactors = FALSE
)

cat("\n--- Running Multi-Variable Limma Models ---\n")
for (g in groupings) {
  cat("Processing Limma for:", g, "... ")
  
  # Extract operational subsets filtering out missing rows per model matrix iteration
  valid_samples <- !is.na(meta_stats[[g]]) & !is.na(meta_stats$Age) & !is.na(meta_stats$Sex)
  meta_sub <- meta_stats[valid_samples, ]
  mtx_sub <- mtx_stats[valid_samples, ]
  
  # Purge ghost attributes from missing samples using droplevels
  meta_sub <- droplevels(meta_sub)
  
  if (length(unique(na.omit(meta_sub[[g]]))) < 2) {
    cat("Skipped (less than 2 groups remain).\n")
    next
  }
  
  # Drop flat/uninformative metrics displaying zero variance within the current patient block
  feature_variances <- apply(mtx_sub, 2, var, na.rm = TRUE)
  mtx_sub <- mtx_sub[, feature_variances > 0, drop = FALSE]
  
  # Transpose clean data to conform with limma structural requirements (Metabolites=Rows)
  mtx_limma_sub <- t(mtx_sub)
  
  # Construct multivariable regression formula adjustments to control for confounding parameters
  if (g == "Age") { 
    formula_str <- "~ Age + Sex" 
  } else if (g == "Sex") { 
    formula_str <- "~ Sex + Age" 
  } else { 
    formula_str <- paste("~", g, "+ Age + Sex") 
  }
  
  design <- model.matrix(as.formula(formula_str), data = meta_sub)
  
  # Map design targets
  target_cols <- grep(g, colnames(design))
  if(length(target_cols) == 0) {
    cat("Skipped (no target columns found).\n")
    next
  }
  
  # Execute linear modeling utilizing generalized generalized least squares tracking
  fit <- lmFit(mtx_limma_sub, design)
  
  # Discard non-estimable parameters displaying structural collinearity/missingness
  estimable_targets <- target_cols[colSums(!is.na(fit$coefficients[, target_cols, drop=FALSE])) > 0]
  if(length(estimable_targets) == 0) {
    cat("Skipped (coefficients not estimable).\n")
    next
  }
  
  # Apply Empirical Bayes moderation to borrow variance across all features, stabilizing small group profiles
  fit <- eBayes(fit)
  
  # Tabulate robust linear modeling metrics with Benjamini-Hochberg (BH) false discovery rate control
  trait_top_table <- topTable(fit, coef = estimable_targets, number = Inf, adjust.method = "BH")
  trait_top_table$Metabolite_ID <- rownames(trait_top_table)
  
  # Index comprehensive output results to master list object
  all_limma_results[[g]] <- trait_top_table
  
  # Quantify statistical tracking thresholds across distinct confidence dimensions
  raw_05 <- sum(trait_top_table$P.Value < 0.05, na.rm = TRUE)
  raw_01 <- sum(trait_top_table$P.Value < 0.01, na.rm = TRUE)
  fdr_05 <- sum(trait_top_table$adj.P.Val < 0.05, na.rm = TRUE)
  
  limma_summary_table <- rbind(limma_summary_table, data.frame(
    Clinical_Trait = g, Features_Raw_P_05 = raw_05, Features_Raw_P_01 = raw_01, Features_Strict_FDR = fdr_05
  ))
  cat("Done.\n")
}

print("BIOMARKER YIELD PER CLINICAL TRAIT (CONTROLLED FOR AGE & SEX):")
print(limma_summary_table)
write.table (limma_summary_table, "~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metabolic_Project/RP_Neg/Limma_summary_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# --- Auto-Export Functional Result Output Tables ---
for (g in groupings) {
  if (!is.null(all_limma_results[[g]])) {
    # paste0 generates clean string connections avoiding file-naming path spaces
    file_name <- paste0("~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metabolic_Project/RP_Neg/Stats_All_Metabolites_", g, ".txt")
    write.table(all_limma_results[[g]], file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

# --- Visualizing Phase II Summaries (Multivariate Yield Plots) ---
# 1. Reshape the summary metrics table into long format for ggplot grouped layout
summary_long <- tidyr::pivot_longer(limma_summary_table, 
                                    cols = starts_with("Features_"), 
                                    names_to = "Strictness", 
                                    values_to = "Count")

# 2. Re-level and rename significance tiers to create a clear logical hierarchy
summary_long$Strictness <- factor(summary_long$Strictness, 
                                  levels = c("Features_Strict_FDR", "Features_Raw_P_01", "Features_Raw_P_05"),
                                  labels = c("FDR < 0.05 (Strict Hits)", "Raw p < 0.01 (Strong Leads)", "Raw p < 0.05 (Exploratory)"))

# 3. CRITICAL: Rank the Clinical Traits on the axis by total significant feature volume
# This ensures your highest-yielding clinical traits sit proudly at the top of the figure
trait_order <- limma_summary_table %>%
  mutate(Total_Signal = Features_Raw_P_05 + Features_Raw_P_01 + Features_Strict_FDR) %>%
  arrange(Total_Signal) %>%
  pull(Clinical_Trait)

summary_long$Clinical_Trait <- factor(summary_long$Clinical_Trait, levels = trait_order)

# 4. Define an ultra-premium editorial color profile
yield_colors <- c(
  "FDR < 0.05 (Strict Hits)"  = "#D90429", # Crimson Red (Draws focus to high-confidence hits)
  "Raw p < 0.01 (Strong Leads)" = "#2B2D42", # Deep Slate Slate Navy (For solid statistical leads)
  "Raw p < 0.05 (Exploratory)"  = "#8D99AE"  # Muted Silver Blue (Soft background layer for broad signals)
)

# 5. Build the Enhanced Grouped Yield Plot
p_yield <- ggplot(summary_long, aes(x = Clinical_Trait, y = Count, fill = Strictness)) +
  # Draw crisp, thin-bordered bars with clean spacing gaps between groups
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           color = "black", width = 0.7, alpha = 0.9) +
  
  # SMART LABELS: Dynamically overlays count values next to active bars
  # (ifelse(Count > 0...) completely suppresses '0' text strings to eliminate baseline noise)
  geom_text(aes(label = ifelse(Count > 0, Count, "")), 
            position = position_dodge(width = 0.8), hjust = -0.25, vjust = 0.4,
            family = "serif", fontface = "bold", size = 3.8, color = "grey10") +
  
  # Implement our custom editorial color scheme
  scale_fill_manual(values = yield_colors, name = "Significance Screening Tier") +
  
  # Flip coordinates for effortless left-to-right reading of clinical categories
  coord_flip() +
  
  # Apply clean minimalist canvas setting
  theme_minimal() +
  labs(
    title = "Differential Regulation Yield Across Patient Phenotypes (RP Negative)",
    subtitle = "Quantified feature discovery rates via multi-variable Limma models (Adjusted for Age & Sex)",
    x = "Evaluated Clinical Phenotype / Trait",
    y = "Number of Differentially Regulated Metabolites (Features Count)"
  ) +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(size = 15, face = "bold", color = "black", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 11, face = "italic", color = "grey40", margin = margin(b = 15)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "top", # Places legend cleanly at top to save horizontal space
    legend.box.spacing = unit(0.3, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank() # Strips horizontal lines to let the trait rows breathe
  ) +
  
  # Expand upper coordinate limits to ensure floating text values never get cut off by the edge
  scale_y_continuous(
    breaks = seq(0, max(summary_long$Count, na.rm = TRUE) + 1, by = 5),
    expand = expansion(mult = c(0, 0.15))
  )

# 6. Render the finalized publication graphic
print(p_yield)

# -------------------------------------------------------------------------
# MODULE 7: TARGETED ADAPTIVE GRAPHICS FOR CHOSEN FOCUS PHENOTYPES
# -------------------------------------------------------------------------
# Global Significance Threshold
p_cutoff <- 0.05

# Define target folder path (space-free naming format)
output_dir <- "~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metabolic_Project/RP_Neg/"

# Custom hyphen/slash-aware text wrapper for long chemical name strings
wrap_chem_names <- function(names, width = 22) {
  sapply(names, function(x) {
    spaced <- gsub("([-–\\/])", "\\1 ", x)
    wrapped <- paste(strwrap(spaced, width = width), collapse = "\n")
    cleaned <- gsub("([-–\\/]) ", "\\1", wrapped)
    cleaned <- gsub("([-–\\/])\n ", "\\1\n", cleaned)
    return(cleaned)
  })
}

cat("\n--- Starting Combined Automated Panel Generation Pipeline (A4 Format) ---\n")

for (g in groupings) {
  cat("Generating composite plate for trait:", g, "... ")
  
  # 1. Verification Check: Ensure the trait has a corresponding limma table
  if (is.null(all_limma_results[[g]])) {
    cat("SKIPPED (No statistical model table found)\n")
    next
  }
  
  # Create a local copy to protect the master list data arrays from mutation
  plot_data <- all_limma_results[[g]]
  
  # =========================================================================
  # STEP A: GENERATE THE LEFT PANEL (MACRO SCREENING PROFILE)
  # =========================================================================
  if ("logFC" %in% colnames(plot_data)) {
    
    # --- OPTION A.1: FOUR-QUADRANT VOLCANO PLOT (FOR BINARY TRAITS) ---
    plot_data$Quadrant <- "Non-significant Down Regulated"
    plot_data$Quadrant[plot_data$P.Value >= p_cutoff & plot_data$logFC > 0] <- "Non-significant Up Regulated"
    plot_data$Quadrant[plot_data$P.Value < p_cutoff & plot_data$logFC < 0] <- "Significantly Down Regulated"
    plot_data$Quadrant[plot_data$P.Value < p_cutoff & plot_data$logFC > 0] <- "Significantly Up Regulated"
    
    plot_data$Quadrant <- factor(plot_data$Quadrant, 
                                 levels = c("Significantly Up Regulated", "Significantly Down Regulated", 
                                            "Non-significant Up Regulated", "Non-significant Down Regulated"))
    
    left_colors <- c(
      "Significantly Up Regulated"     = "#E69F00", 
      "Significantly Down Regulated"   = "#2563EB", 
      "Non-significant Up Regulated"   = "#F3E1D3", 
      "Non-significant Down Regulated" = "#DBE2EF"  
    )
    
    label_hits <- plot_data %>% 
      filter(P.Value < p_cutoff) %>% 
      arrange(P.Value) %>% 
      head(10)
    
    p_left <- ggplot(plot_data, aes(x = logFC, y = -log10(P.Value))) +
      geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey60", linewidth = 0.4) +
      geom_vline(xintercept = 0, linetype = "solid", color = "grey40", linewidth = 0.4) +
      geom_point(aes(color = Quadrant), alpha = 0.85, size = 2.4, stroke = 0.2) +
      scale_color_manual(values = left_colors, name = "Metabolic Regulation") +
      geom_text_repel(
        data = label_hits, aes(label = Metabolite_ID),
        size = 3.6, fontface = "bold", color = "black", family = "serif",
        box.padding = 0.5, point.padding = 0.3, max.overlaps = 20, force = 3,
        segment.size = 0.35, segment.color = "grey30", arrow = arrow(length = unit(0.015, "npc"))
      ) +
      labs(
        title = paste("Panel A: Volcano Trajectory Profiling :", g), 
        subtitle = paste("Features partitioned dynamically at raw p <", p_cutoff, "(Controlled for Age & Sex)"),
        x = expression(bold("Directional Effect Size (Log"[2]*" Fold Change)")),
        y = expression(bold("-Log"[10]*" Raw P-value"))
      )
    
    file_prefix <- "Composite_Plate_Volcano_"
    
  } else if ("F" %in% colnames(plot_data)) {
    
    # --- OPTION A.2: MULTI-CLASS SIGNIFICANCE TRACK F-PLOT (FOR COMPLEX TRAITS) ---
    plot_data$Significance_Tier <- "Not Significant"
    plot_data$Significance_Tier[plot_data$P.Value < 0.05] <- "Significant (p < 0.05)"
    plot_data$Significance_Tier[plot_data$P.Value < 0.01] <- "Highly Significant (p < 0.01)"
    
    plot_data$Significance_Tier <- factor(plot_data$Significance_Tier, 
                                          levels = c("Highly Significant (p < 0.01)", "Significant (p < 0.05)", "Not Significant"))
    
    left_colors <- c(
      "Highly Significant (p < 0.01)" = "#6B21A8", 
      "Significant (p < 0.05)"        = "#B06AB3", 
      "Not Significant"               = "#E5E7EB"  
    )
    
    p_left <- ggplot(plot_data, aes(x = F, y = -log10(P.Value))) +
      geom_hline(yintercept = c(-log10(0.05), -log10(0.01)), linetype = "dashed", color = "grey60", alpha = 0.5) +
      geom_point(aes(color = Significance_Tier), alpha = 0.8, size = 2.5) +
      scale_color_manual(values = left_colors, name = "Statistical Weight") +
      geom_text_repel(
        data = head(plot_data[order(plot_data$P.Value), ], 10), aes(label = Metabolite_ID),
        size = 3.6, fontface = "bold", color = "black", family = "serif",
        box.padding = 0.5, point.padding = 0.3, max.overlaps = 15, force = 2,
        segment.size = 0.35, segment.color = "grey40", arrow = arrow(length = unit(0.015, "npc"))
      ) +
      labs(
        title = paste("Panel A: Multi-Group Variance Profile :", g), 
        subtitle = "Identifies features shifting significantly across any group level (Controlled for Age & Sex)",
        x = "Overall Group Variance Weight (Limma F-Statistic)",
        y = expression(bold("-Log"[10]*" Raw P-value"))
      )
    
    file_prefix <- "Composite_Plate_MultiGroup_"
  }
  
  p_left <- p_left + 
    theme_classic() + 
    theme(
      text = element_text(family = "serif"),
      plot.title = element_text(face = "bold", size = 13, color = "grey10", margin = margin(b = 4)),
      plot.subtitle = element_text(size = 9, face = "italic", color = "grey40", margin = margin(b = 15)),
      axis.title.x = element_text(size = 11, color = "black", margin = margin(t = 10)),
      axis.title.y = element_text(size = 11, color = "black", margin = margin(r = 10)),
      axis.text = element_text(size = 10, color = "grey10"),
      axis.line = element_line(linewidth = 0.5, color = "grey20"), 
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.position = "bottom"
    )
  
  # =========================================================================
  # STEP B: GENERATE THE RIGHT PANEL (BEAUTIFIED TOP 4 TRAJECTORY BOXPLOT GRID)
  # =========================================================================
  top_features <- plot_data %>% arrange(P.Value) %>% head(4) %>% pull(Metabolite_ID)
  
  # MODIFICATION 1: Allow "Short Lasting Arthritis" controls to remain in the plot ONLY for Disease_Type
  if (g %in% c("Disease_activity", "Disease_Type", "Disease_outcome", "TMJ_arthritis", "Cervical_arthritis")) {
    if (g == "Disease_Type") {
      cohort_filter <- rep(TRUE, nrow(meta_stats)) # Include controls for baseline comparison
    } else {
      cohort_filter <- meta_stats$Group != "Short Lasting Arthritis"
    }
  } else {
    cohort_filter <- rep(TRUE, nrow(meta_stats)) 
  }
  
  valid_rows <- cohort_filter & !is.na(meta_stats[[g]])
  
  grid_df <- data.frame(mtx_stats[valid_rows, top_features, drop = FALSE], check.names = FALSE)
  grid_df$Group_Level <- as.factor(meta_stats[[g]][valid_rows])
  
  grid_long <- tidyr::pivot_longer(grid_df, cols = all_of(top_features), 
                                   names_to = "Metabolite", values_to = "Intensity")
  
  # MODIFICATION 2: Map factor level "0" to "Short Lasting" for Disease_Type
  if (g == "Disease_Type") {
    grid_long$Group_Level <- factor(grid_long$Group_Level, levels = c("0", "1", "2", "3"), labels = c("SLA", "ERA", "Oligo", "Poly"))
  } else if (g == "Disease_activity") {
    grid_long$Group_Level <- factor(grid_long$Group_Level, levels = c("1", "2", "3"), labels = c("Low", "Moderate", "High"))
  } else if (g == "Sex") {
    grid_long$Group_Level <- factor(grid_long$Group_Level, levels = c("1", "2"), labels = c("Female", "Male"))
  } else if (g == "Age") {
    grid_long$Group_Level <- factor(grid_long$Group_Level, levels = c("0", "1"), labels = c("1-4.9 yrs", "4.9+ yrs"))
  } else {
    grid_long$Group_Level <- factor(grid_long$Group_Level, levels = c("0", "1"), labels = c("Negative/No", "Positive/Yes"))
  }
  
  # MODIFICATION 3: Define custom 4-group color narrative (Neutral Gray for baseline Control)
  num_groups <- length(levels(grid_long$Group_Level))
  if (num_groups == 2) {
    box_palette <- c("#1A5F7A", "#57C5B6") 
  } else if (num_groups == 3) {
    box_palette <- c("#2A2F4F", "#917FB3", "#E5BEEC") 
  } else if (num_groups == 4) {
    box_palette <- c("#8D99AE", "#2A2F4F", "#917FB3", "#E5BEEC") # Gray-blue Control, Amethyst tones for disease subtypes
  } else {
    box_palette <- c("#0F2C59", "#DAC0A3", "#E48585", "#7091F5")
  }
  
  p_right <- ggplot(grid_long, aes(x = Group_Level, y = Intensity, fill = Group_Level)) +
    geom_boxplot(alpha = 0.50, outlier.shape = NA, width = 0.40, 
                 color = "grey20", linewidth = 0.45, median.linewidth = 1.1) +
    
    geom_jitter(aes(fill = Group_Level), shape = 21, color = "white", 
                width = 0.10, size = 1.8, stroke = 0.4, alpha = 0.8) +
    
    facet_wrap(~Metabolite, scales = "free_y", nrow = 2, 
               labeller = as_labeller(wrap_chem_names)) +
    
    scale_fill_manual(values = box_palette) +
    
    theme_minimal() + 
    labs(
      title = "Panel B: Top 4 Trajectory Resolution Grid",
      subtitle = "True demographic-adjusted distribution dynamics across subcohort levels",
      x = paste("Patient Category:", g),
      y = "Normalized & Scaled Metabolic Intensity"
    ) +
    theme(
      text = element_text(family = "serif"),
      plot.title = element_text(face = "bold", size = 13, color = "grey10", margin = margin(b = 4)),
      plot.subtitle = element_text(size = 9, face = "italic", color = "grey40", margin = margin(b = 15)),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 9, color = "black", hjust = 0.0, margin = margin(b = 5), lineheight = 0.85),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3, color = "grey90"),
      axis.title.x = element_text(size = 11, color = "black", margin = margin(t = 10)),
      axis.title.y = element_text(size = 11, color = "black", margin = margin(r = 10)),
      axis.text.x = element_text(size = 10, face = "bold", color = "grey20"),
      axis.text.y = element_text(size = 9, color = "grey40"),
      legend.position = "none"
    )
  
  # =========================================================================
  # STEP C: PANELS ASSEMBLY VIA WRAP_PLOTS
  # =========================================================================
  p_composite <- patchwork::wrap_plots(p_left, p_right) + patchwork::plot_layout(widths = c(1, 1.25))
  
  full_save_path <- paste0(output_dir, file_prefix, g, ".pdf")
  
  ggsave(
    filename = full_save_path,
    plot = p_composite,
    width = 297,   
    height = 210,  
    units = "mm",
    dpi = 300      
  )
  
  cat("SAVED COMPOSITE PLATE PANELS SUCCESSFULLY.\n")
}

cat("\n--- Pipeline Completed. All integrated A4 plates exported successfully ---\n")
# =========================================================================
# END OF PIPELINE
# ========================================================================
