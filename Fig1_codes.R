# ---- Load Required Libraries ----
library(multcomp)
library(car)
library(ggplot2)       
library(gridExtra)     
library(ggrepel)       
library(ggpubr)      
library(stats) 
library(biostat)

################################################################
# ---- Load and Label Datasets for Analysis ----
load("Fig1.RData")

# Plant biometric measurements (e.g., shoot weight, root length & height)
Data_all      

# Elemental concentrations in root tissue via ICP-OES
ICP_Data_all  

# Soil biochemical traits (e.g., SOM, TC, TN, soil moisture)
Soildata      


#Figure 1B
# Ensure consistent ordering of treatment levels for plotting and contrast coding
Data_all$Treatment <- factor(Data_all$Treatment,
                             levels = c("Control", "EBPR", "Vermont",
                                        "AMF", "EBPR_AMF", "Vermont_AMF"))

# Apply Z-score transformation to all numeric columns
numeric_cols <- sapply(Data_all, is.numeric)
Data_all_z <- as.data.frame(scale(Data_all[, numeric_cols]))
Data_all_z$Treatment<-Data_all$Treatment


# ---- Define Trait List and Output Containers ----
traits <- c("Stalk", "Adventitious", "Shootheight", "Rootlength", "AMFcolonization", "Shootweight")
results <- list()

# ---- Loop Through Each Trait ----
for (trait in traits) {
  cat("\n===============================\n")
  cat(paste("🌿 TRAIT:", trait, "\n"))
  cat("===============================\n")
  
  # Apply post hoc tests
  cat("→ Pairwise Wilcoxon with BH adjustment and CLD:\n")
  
  # Apply pairwise Wilcoxon test with BH correction
  wilcox_res <- pairwise.wilcox.test(Data_all_z[[trait]], 
                                     Data_all_z$Treatment, p.adjust.method = "BH")
  
  # Generate compact letter display
  cld_letters <- biostat::make_cld(wilcox_res)
  
  # Summarize results
  posthoc_summary <- list(kruskal = wilcox_res, cld = cld_letters)
  
  # 5. Save results
  results[[trait]] <- list(posthoc = posthoc_summary)
}
results

results$Adventitious
#Create figures

# ---- Custom Color Palettes for Treatment Groups ----
fill_colors <- c("seagreen3", "cornflowerblue", "orange2", "coral", "royalblue4", "coral4")
color_map   <- fill_colors  # Same for both fill and outline

# ---- Traits to Plot with Corresponding Y-Axis Labels ----
traits <- c("Shootheight", "Shootweight", "Rootlength",
            "Stalk", "Adventitious", "AMFcolonization")
ylabs <- c("Shoot height (cm)",
           "Shoot dry mass (g)",
           "Root length (cm)",
           "Stalk diameter (cm)",
           "Number of adventitious roots",
           "AMF colonization rate (%)")

# ---- Create ggplot Objects in a Loop ----
ggplot_list <- list()

for (i in seq_along(traits)) {
  trait <- traits[i]
  ylabel <- ylabs[i]
  
  # Dynamically construct ggplot for each trait
  p <- ggplot(data = Data_all, aes_string(x = "Treatment", y = trait,
                                          color = "Treatment", fill = "Treatment")) +
    geom_boxplot(aes_string(x = "Treatment", y = trait,
                            fill = "Treatment", color = "Treatment"),
                 outlier.size = 0.4, width = 0.8, alpha = 0.5,
                 position = position_dodge(width = 0.9)) +
    geom_point(aes_string(x = "Treatment", y = trait,
                          fill = "Treatment", color = "Treatment"),
               position = position_jitterdodge(dodge.width = 0.6)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 12),
          legend.position = "none",
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          aspect.ratio = 0.5) +
    theme_classic() +
    labs(y = ylabel) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = color_map) +
    scale_y_continuous(expand = c(0.1, 0, 0.1, 0))
  
  # Store the plot object
  ggplot_list[[trait]] <- p
}

# ---- Arrange Plots in Custom Grid (3x2 Layout) ----
# Layout matrix defines positioning: each number refers to a plot in order
layout_matrix <- rbind(c(1,1,2,2),
                       c(3,3,4,4),
                       c(5,5,6,6))

# Extract plots in order of traits
Plantggplot <- grid.arrange(grobs = ggplot_list[traits],
                            layout_matrix = layout_matrix,
                            respect = TRUE)

#Figure 1C

# Reading a space-separated file without a header
Soildata$Treatment<-factor(Soildata$Treatment, levels = c("Control","EBPR", "Vermont", 
                                                 "AMF","EBPR_AMF", "Vermont_AMF" ))
# Apply Z-score transformation to all numeric columns
numeric_cols <- sapply(Soildata, is.numeric)
Soildata_z <- as.data.frame(scale(Soildata[, numeric_cols]))
Soildata_z$Treatment<-Soildata$Treatment

# ---- Define Target Elements and Output Storage ----
Soildatalist <- c("SOM", "TC", "EGRSP", "TGRSP", "TN", "Moisture")
results_list_soildata <- list()

# Loop through each trait
for (dl in Soildatalist) {
  cat("\n--- Analyzing:", dl, "---\n")
  
  #  Step 1: Compute mean ± SE per treatment
  stats_df <- Soildata %>%
    group_by(Treatment) %>%
    dplyr::summarise(
      mean = mean(.data[[dl]], na.rm = TRUE),
      se = sd(.data[[dl]], na.rm = TRUE) / sqrt(dplyr::n())
    )
  # Step 2: Fit model on z-transformed trait values
  trait_formula <- as.formula(paste(dl, "~ Treatment"))
  lm_model <- lm(trait_formula, data = Soildata_z)
  
  # Step 3: Tukey post-hoc comparison with BH-adjusted p-values
  tukey_test <- glht(lm_model, linfct = mcp(Treatment = "Tukey"))
  tukey_summary <- summary(tukey_test, test = adjusted("BH"))
  
  # Step 4: Generate compact letter display (CLD)
  cld_letters <- cld(tukey_test, test = adjusted("BH"))
  cld_df <- data.frame(
    Treatment = names(cld_letters$mcletters$Letters),
    CLD = cld_letters$mcletters$Letters
  )
  
  # Step 5: Merge CLD with treatment-wise mean ± SE
  summary_df <- dplyr::left_join(stats_df, cld_df, by = "Treatment")
  
  # Step 6: Save annotated results
  results_list_soildata[[dl]] <- list(
    tukey_summary = tukey_summary,
    cld_letters = cld_letters,
    summary_table = summary_df
  )
}


#Figure 1D
ICP_Data_all$Treatment<-factor(ICP_Data_all$Treatment, levels = c("Control","EBPR", "Vermont", 
                                                          "AMF","EBPR_AMF", "Vermont_AMF" ))

# Apply Z-score transformation to all numeric columns
numeric_cols <- sapply(ICP_Data_all, is.numeric)
ICP_Data_all_z <- as.data.frame(scale(ICP_Data_all[, numeric_cols]))
ICP_Data_all_z$Treatment<-ICP_Data_all$Treatment

# ---- Define Target Elements and Output Storage ----
elements <- c("Ca", "K", "Mg", "P", "Fe", "S", "Al", "Zn")
results_list <- list()

# ---- Loop Through Each Trait ----
for (el in elements) {
  cat("\n--- Analyzing:", el, "---\n")
  
  # Apply post hoc tests
  cat("→ Pairwise Wilcoxon with BH adjustment and CLD:\n")
  
  # Apply pairwise Wilcoxon test with BH correction
  wilcox_res <- pairwise.wilcox.test(ICP_Data_all_z[[el]],
                                     ICP_Data_all_z$Treatment, p.adjust.method = "BH")
  
  # Generate compact letter display
  cld_letters <- biostat::make_cld(wilcox_res)
  
  # Summarize results
  posthoc_summary <- list(kruskal = wilcox_res, cld = cld_letters)
  
  # Save results
  results_list[[el]] <- list(
    posthoc_summary = posthoc_summary,
    cld_letters = cld_letters
  )
}

minerals <- c("Al", "Ca", "Fe", "K", "Mg", "P", "S", "Zn")  # Minerals to visualize
y_labels <- pCay_labels <- paste(minerals, "(mg/g)")                      # Y-axis labels
fill_colors <- c("seagreen3", "cornflowerblue", "orange2", "coral",
                 "royalblue4", "coral4")                   # Consistent treatment palette

# ---- Create Plot List via Loop ----
mineral_plots <- list()

for (i in seq_along(minerals)) {
  mineral <- minerals[i]
  ylabel <- y_labels[i]
  
  p <- ggplot(ICP_Data_all, aes_string(x = "Treatment", y = mineral,
                                            color = "Treatment", fill = "Treatment")) +
    geom_boxplot(outlier.size = 0.4, width = 0.8, alpha = 0.5,
                 position = position_dodge(width = 0.9)) +
    geom_point(position = position_jitterdodge(dodge.width = 0.9)) +
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Adds border
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 12),
      legend.position = "none",
      aspect.ratio = 4
    ) +
    labs(y = ylabel) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors); p
  
  mineral_plots[[mineral]] <- p
}

# ---- Arrange All Plots in a 4×2 Layout ----
layout_matrix <- rbind(c(1, 1, 2, 2, 3, 3, 4, 4),
                       c(1, 1, 2, 2, 3, 3, 4, 4),
                       c(5, 5, 6, 6, 7, 7, 8, 8),
                       c(5, 5, 6, 6, 7, 7, 8, 8))
ICPggplot <- grid.arrange(grobs = mineral_plots[minerals],
                          layout_matrix = layout_matrix,
                          respect = TRUE)

