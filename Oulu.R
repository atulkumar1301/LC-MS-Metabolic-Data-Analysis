library(data.table)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(car)
library(emmeans)
library(broom)
library(visreg)
library(ggsignif)
library(patchwork)
library(stringr)
library(UpSetR)
library(dplyr)
library(tidyr)
library(cowplot)
library(gridExtra)
options(scipen = 999)

df <- fread ("~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metabolic_Project/RP_Pos/3_Clean_Data.txt")

mtx <- df [, 4:68]
unique_names <- make.unique(as.character(df[[1]]))
rownames(mtx) <- unique_names

# Median Normalization (Vertical)
sample_medians <- apply(mtx, 2, median, na.rm = TRUE)
mtx_norm <- sweep(mtx, 2, sample_medians, "/") * median(sample_medians)

# Log Transformation (To handle the skew)
mtx_log <- log2(mtx_norm + 1) # Adding 1 avoids log(0) issues

# Pareto Scaling (Horizontal)
# This makes small-intensity peaks visible next to big ones
mtx_final <- t(apply(mtx_log, 1, function(x) (x - mean(x)) / sqrt(sd(x))))

# Load metadata
# It should have columns: SampleID and Group
metadata <- fread("~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metadata_pos_RP.txt")

# Run PCA
pca_res <- prcomp(t(mtx_final), center = FALSE, scale. = FALSE)

# Calculate variance explained for the axes
var_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

# Create a data frame for plotting
pca_df <- data.frame(SampleID = rownames(pca_res$x),
                     PC1 = pca_res$x[,1],
                     PC2 = pca_res$x[,2],
                     Group = metadata$Group)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA: Metabolic Profiling (Positive RP)",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  scale_color_brewer(palette = "Set1") + # Good for 4-5 groups
  stat_ellipse(level = 0.95) # Draws 95% confidence circles around groups
p <- p + scale_x_continuous(breaks = round(seq(-50, 40, by = 5),1))
p <- p + scale_y_continuous(breaks = round (seq (-25, 20, by = 5), 1))
p <- p +
  theme(plot.title = element_text(family = "serif", size=18, face = "bold"),
        axis.title.x = element_text(family = "serif", size=16),
        axis.title.y = element_text(family = "serif", size=16),
        axis.text.x = element_text(family = "serif", size=14),
        axis.text.y = element_text(family = "serif", size=14),
        legend.title = element_text(family = "serif", size=16),
        legend.text = element_text(family = "serif", size=16),
        panel.background = element_blank())

p

# We only want to compare the 4 biological groups so remove QC data
norm_data <- as.data.frame(t(mtx_final))
is_qc <- metadata$Group == "QC" 
mtx_stats <- norm_data[!is_qc, ]
meta_stats <- metadata[!is_qc, ]

results_list <- list()

for (metab in colnames(mtx_stats)){
  # Build local dataframe
  df_1 <- data.frame(
    val = mtx_stats[[metab]],
    Group = meta_stats$Group,
    Age = meta_stats$Age,
    Sex = meta_stats$Gender)
  
  # Fit ANCOVA (Type II)
  fit <- lm(val ~ Group + Age + Sex, data = df_1)
  anova_res <- Anova(fit, type = "II")
  
  # Calculate Partial Eta Squared (Effect Size)
  ss_group <- anova_res["Group", "Sum Sq"]
  ss_resid <- anova_res["Residuals", "Sum Sq"]
  eta_sq   <- ss_group / (ss_group + ss_resid)
  
  # Pairwise Tukey Comparisons
  emm <- emmeans(fit, specs = "Group")
  tukey_df <- as.data.frame(pairs(emm, adjust = "tukey"))
  
  # Create a named vector for this metabolite's stats
  stats_vec <- c(
    Metabolite      = metab,
    Group_F         = anova_res["Group", "F value"],
    Group_P         = anova_res["Group", "Pr(>F)"],
    Effect_Size_Eta2 = eta_sq,
    Age_P           = anova_res["Age", "Pr(>F)"],
    Sex_P           = anova_res["Sex", "Pr(>F)"]
  )
  
  # Dynamically add Tukey results based on contrast names
  for(i in 1:nrow(tukey_df)) {
    contrast_name <- gsub(" ", "", tukey_df$contrast[i])
    stats_vec[paste0(contrast_name, "_p_adj")] <- tukey_df$p.value[i]
    stats_vec[paste0(contrast_name, "_diff")]  <- tukey_df$estimate[i]
  }
  
  results_list[[metab]] <- stats_vec
}

# Combine into final table
TABLE <- as.data.frame(do.call(rbind, results_list))

# Convert numeric columns (they become characters during rbind)
numeric_cols <- colnames(TABLE)[colnames(TABLE) != "Metabolite"]
TABLE[numeric_cols] <- lapply(TABLE[numeric_cols], as.numeric)

write.table (TABLE, "~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metabolic_Project/RP_Pos/4_ANOVA_Outcome.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)  


# Function to generate a grid of plots based on a filtered results table
plot_metabolites <- function(sig_table, original_mtx, metadata, ncol = 4) {
  metabolites <- sig_table$Metabolite
  if(length(metabolites) == 0) return(message("No significant metabolites found."))
  
  plot_list <- list()
  
  for (metabo in metabolites) {
    # Prepare data for current metabolite
    df_plot <- data.frame(
      val = original_mtx[[metabo]],
      Group = metadata$Group,
      Age = metadata$Age,
      Gender = metadata$Gender
    )
    
    # Fit ANCOVA and extract adjusted residuals (visreg handles Age/Sex adjustment)
    fit <- lm(val ~ Group + Age + Gender, data = df_plot)
    v <- visreg(fit, "Group", plot = FALSE)
    adj_df <- v$res
    
    # Generate Plot
    p <- ggplot(adj_df, aes(x = Group, y = visregRes, fill = Group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.1, alpha = 0.4, size = 0.8) +
      theme_classic(base_family = "serif") +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = str_wrap(metabo, width = 15), y = "Adj. Intensity", x = "Groups")
    p <- p +
      theme(plot.title = element_text(size=8, face = "bold"),
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12),
            axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=10),
            legend.title = element_text(size=12),
            legend.text = element_text(size=12),
            plot.margin = margin(t = 5, r = 10, b = 5, l = 10),
            panel.background = element_blank())
    
    plot_list[[metabo]] <- p
  }
  
  #return(wrap_plots(plot_list, ncol = ncol))
  return (plot_list)
}


# Define comparisons: Column Name = File Suffix
comparisons <- list(
  "Group_P"                            = list(sfx = "ANCOVA",        nc = 3),
  "Control-JIA_p_adj"                  = list(sfx = "Control_JIA",   nc = 3),
  "Control-Oligoarthritis_p_adj"       = list(sfx = "Control_Oligo",  nc = 2),
  "Control-Polyarthritis_p_adj"        = list(sfx = "Control_Poly",   nc = 1),
  "JIA-Oligoarthritis_p_adj"           = list(sfx = "JIA_Oligo",     nc = 1),
  "JIA-Polyarthritis_p_adj"            = list(sfx = "JIA_Poly",      nc = 3),
  "Oligoarthritis-Polyarthritis_p_adj" = list(sfx = "Oligo_Poly",     nc = 2)
)

base_path <- "~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metabolic_Project/RP_Pos/Box_Plot_Significant_Metabolite_"

for (col in names(comparisons)) {
  sub_table <- TABLE[TABLE[[col]] < 0.05, ]
  
  if (nrow(sub_table) > 0) {
    plots <- plot_metabolites(sub_table, mtx_stats, meta_stats)
    
    # Extract specific config for this comparison
    suffix <- comparisons[[col]]$sfx
    n_col  <- comparisons[[col]]$nc
    
    pdf(paste0(base_path, suffix, ".pdf"), width = 11.69, height = 8.27)
    
    # Using grid.draw logic or marrangeGrob directly
    grid_plots <- print(marrangeGrob(plots, nrow = 3, ncol = n_col, top = NULL))
    dev.off()
  } else {
    message(paste("No significant metabolites for", col))
  }
}

# Heat Map
# Select top metabolites from your TABLE
top_metabs <- na.omit(TABLE[order(TABLE$Group_P), "Metabolite"])[1:20]

# Extract and force to numeric matrix
hm_mtx <- t(mtx_stats[, top_metabs])
hm_mtx <- matrix(as.numeric(hm_mtx), nrow=nrow(hm_mtx), dimnames=list(rownames(hm_mtx), colnames(hm_mtx)))

# THE FIX: Remove ANY row with NA or Zero Variance
# This is the 'Nuclear Option' to ensure pheatmap has finite numbers
clean_rows <- apply(hm_mtx, 1, function(x) {
  finite <- all(is.finite(x))        # No NAs or Infs
  variation <- sd(x, na.rm=TRUE) > 0 # Must have some change
  return(finite & variation)
})

hm_mtx_final <- hm_mtx[clean_rows, ]


# Check if we have anything left
if(nrow(hm_mtx_final) == 0) {
  stop("Error: No metabolites passed the finite/variance check. Check your mtx_stats values.")
}

# Ensure metadata matches
anno_col <- data.frame(
  Group = meta_stats$Group,
  Age = meta_stats$Age,
  Sex = meta_stats$Gender
)
rownames(anno_col) <- colnames(hm_mtx_final)

# This replaces spaces with newlines every 15 chars
rownames(hm_mtx_final) <- str_wrap(rownames(hm_mtx_final), width = 15)

pheatmap(hm_mtx_final, 
         annotation_col = anno_col, 
         scale = "row",                # Now safe because variance > 0
         clustering_distance_rows = "correlation",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         show_colnames = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 5,             # Small font for wrapped text
         lineheight = 0.8,            # Tighten the space between wrapped lines
         #cellwidth = 20,              # Manual cell size prevents shrinking
         #cellheight = 15,             # Ensures rows are tall enough for wrapped text
         main = "Heatmap Top 20 Metabolites (RP Positive)")

# Upset-Plot

Anova_sig <- TABLE [TABLE$Group_P < 0.05, "Metabolite"]
Control_JIA_sig <- TABLE [TABLE$`Control-JIA_p_adj` < 0.05, "Metabolite"]
Control_Oligo_sig <- TABLE [TABLE$`Control-Oligoarthritis_p_adj` < 0.05, "Metabolite"]
Control_Poly_sig <- TABLE [TABLE$`Control-Polyarthritis_p_adj` < 0.05, "Metabolite"]
Jia_Oligo_sig <- TABLE [TABLE$`JIA-Oligoarthritis_p_adj` < 0.05, "Metabolite"]
Jia_Poly_sig <- TABLE [TABLE$`JIA-Polyarthritis_p_adj` < 0.05, "Metabolite"]
Oligo_Poly_sig <- TABLE [TABLE$`Oligoarthritis-Polyarthritis_p_adj` < 0.05, "Metabolite"]

venn_list <- list(ANCOVA = na.omit(Anova_sig), 
                  Control_JIA = na.omit(Control_JIA_sig), 
                  Control_Oligo = na.omit(Control_Oligo_sig), 
                  Control_Poly = na.omit(Control_Poly_sig),
                  JIA_Oligo = na.omit(Jia_Oligo_sig),
                  JIA_Poly = na.omit(Jia_Poly_sig),
                  Oligo_Poly = na.omit(Oligo_Poly_sig))

upset_input <- fromList(venn_list)

upset(upset_input, 
      nsets = 7,                   # Number of groups (Control, JIA, Oligo, Poly)
      order.by = "freq",           # Shows largest intersections first
      decreasing = TRUE,           
      mb.ratio = c(0.6, 0.4),      # Balance between bar chart and dot matrix
      number.angles = 0, 
      text.scale = c(1.3, 1.3, 1, 1, 1.5, 1.2), # Adjust font sizes for labels/axes
      sets.bar.color = "#56B4E9")

# Clean the list: Force every element to be a simple character vector
# This strips away problematic hidden attributes or list structures
clean_list <- lapply(venn_list, function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  return(as.character(unlist(x)))
})

# Filter out groups that have 0 metabolites
final_list <- clean_list[sapply(clean_list, length) > 0]

# Handle the case where the whole list is empty
if (length(final_list) == 0) {
  stop("The list is empty. No metabolites met your significance criteria.")
}

# Manually build the table (Instead of using stack())
venn_df <- data.frame(
  Metabolite = unname(unlist(final_list)),
  Group = rep(names(final_list), sapply(final_list, length)),
  stringsAsFactors = FALSE
)

write.table (venn_df, "~/OneDrive - University of Eastern Finland/Projects/Oulu_Project/Metabolic_Project/RP_Pos/5_Significant_Metabolites.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

graphics.off()

## Volcano Plot

# Define the reusable function
create_volcano <- function(data, diff_col, p_col, group_name1, group_name2, title_text, x_limits = c(-3, 3), y_limits = c(0, 4), x_step = 0.5, y_step = 0.5) {
  
  cbbPalette <- c("#D55E00", "#56B4E9", "#009E73", "#E69F00", "#0072B2", "#F0E442", "#CC79A7", "#000000")
  
  # Prepare data
  df <- data %>%
    select(Metabolite, diff = !!sym(diff_col), p_adj = !!sym(p_col)) %>%
    mutate(
      Regulation = case_when(
        diff > 0 & p_adj < 0.05 ~ paste("Significant Increased", group_name1),
        diff < 0 & p_adj < 0.05 ~ paste("Significant Increased", group_name2),
        diff > 0 & p_adj > 0.05 ~ paste("Non-Significant Increased", group_name1),
        TRUE ~ paste("Non-Significant Increased", group_name2)
      ),
      Label = case_when(
        p_adj < 0.05 ~ as.character(Metabolite),
        TRUE ~ NA_character_
      )
    )
  
  # Build Plot
  ggplot(df, aes(x = diff, y = -log10(p_adj), col = Regulation, label = Label)) +
    geom_point() +
    geom_text_repel(size = 2,                 # Smaller font size (default is usually ~4)
                    max.overlaps = Inf, 
                    force = 50,                 # Higher force pushes labels away from points
                    box.padding = 0.5,          # Small padding around each label
                    point.padding = 0.2,        # Space between label and the point
                    segment.size = 0.2,         # Thinner connector lines
                    min.segment.length = 0,     # Always show lines, even for close labels
                    family = "serif",
                    show.legend = FALSE) +
    geom_hline(aes(yintercept = -log10(0.05), linetype = "p-value 0.05"), col = "black") +
    scale_linetype_manual(name = "p-value cut off", values = 2) +
    scale_color_manual(values = cbbPalette) +
    scale_x_continuous(limits = x_limits, breaks = seq(x_limits[1], x_limits[2], by = x_step)) +
    scale_y_continuous(limits = y_limits, breaks = seq(y_limits[1], y_limits[2], by = y_step)) +
    labs(
      title = title_text,
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~(P[adj])),
      color = "Regulation"
    ) +
    theme_light() +
    theme(
      legend.position = "right", # Removed individual legends for a cleaner grid
      plot.title = element_text(family = "serif", size = 12, face = "bold", hjust = 0.5),
      axis.title = element_text(family = "serif", size = 10),
      axis.text = element_text(family = "serif", size = 10),
      legend.title = element_text(family = "serif", size=10),
      legend.text = element_text(family = "serif", size=10),
      panel.background = element_blank()
    )
}

# Generate the plots using the function
p_C_jia <- create_volcano(TABLE, "Control-JIA_diff", "Control-JIA_p_adj", "Control", "JIA", "Control vs JIA", x_limits = c(-3, 3), y_limits = c(0, 3), x_step = 1.0)
p_C_oligo <- create_volcano(TABLE, "Control-Oligoarthritis_diff", "Control-Oligoarthritis_p_adj", "Control", "Oligoarthritis", "Control vs Oligoarthitis", x_limits = c(-2, 2), y_limits = c(0, 2), x_step = 1.0)
p_C_Poly <- create_volcano(TABLE, "Control-Polyarthritis_diff", "Control-Polyarthritis_p_adj", "Control", "Polyarthritis", "Control vs Polyarthritis", x_limits = c(-2, 2), y_limits = c(0, 2), x_step = 1.0)
p_J_Oligo <- create_volcano(TABLE, "JIA-Oligoarthritis_diff", "JIA-Oligoarthritis_p_adj", "JIA", "Oligoarthritis", "JIA vs Oligoarthritis", x_limits = c(-3, 3), y_limits = c(0, 4), x_step = 1.0, y_step = 1)
p_J_Poly <- create_volcano(TABLE, "JIA-Polyarthritis_diff", "JIA-Polyarthritis_p_adj", "JIA", "Polyarthritis", "JIA vs Polyarthritis", x_limits = c(-3, 3), y_limits = c(0, 3), x_step = 1.0, y_step = 1)
p_O_Poly <- create_volcano(TABLE, "Oligoarthritis-Polyarthritis_diff", "Oligoarthritis-Polyarthritis_p_adj", "Oligoarthritis", "Polyarthritis", "Oligoarthritis vs Polyarthritis", x_limits = c(-1, 1), y_limits = c(0, 1.5))

# Combine into a grid
# Note: I'm adding 'rel_widths' to give space if you want a legend on the side

plot_grid(p_C_jia, p_C_oligo, p_C_Poly, p_J_Oligo, p_J_Poly, p_O_Poly, labels = "AUTO", ncol = 2)


