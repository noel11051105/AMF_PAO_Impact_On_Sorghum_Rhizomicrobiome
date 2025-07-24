

# Load necessary libraries
library(Rfast)
library(ALDEx2)
library(metagMisc)
library(RColorBrewer)
library(ggplot2)
library(phyloseq)
library(vegan)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(microeco)
library(pheatmap)
library(ggtree)
library(ggpicrust2)
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(viridis)

# ---- Loaded Analysis Datasets ----
load("Fig4.RData")

Aldex_glom_genus   # Genus-level ASV data formatted for differential abundance testing with ALDEx2
RA_data            # Relative abundance table of ASVs for taxonomic composition analysis
dataset            # Microeco-compatible microbial abundance dataset for functional or ecological modeling
Vegan_diversity    # Alpha diversity metrics (e.g., Shannon, Simpson) computed using the 'vegan' package
abundance          # Functional KO abundance matrix for pathway-level interpretation
metadata           # Sample metadata including experimental treatments and environmental parameters
temp               # Elemental concentration profiles acquired via ICP-OES
Picrust_result     # PICRUSt2 output with predicted metagenomic functional pathways
KO.tsv             # KO abundance data


# Figure 4A

# ---- Step 1: Filter taxa by ≥70% prevalence across each location group ----
phyloseq_list = list()
for (i in unique(sample_data(Aldex_glom_genus)$Location)) {
  print(i)
  Individualdata <- subset_samples(Aldex_glom_genus, Location %in% i)
  tmp_data <- phyloseq_filter_prevalence(Individualdata, prev.trh = 0.70)
  phyloseq_list[[i]] = tmp_data
}

# ---- Step 2: Merge all location-specific filtered objects into one phyloseq ----
glom_genus2 <- merge_phyloseq(phyloseq_list[[1]], phyloseq_list[[2]])
for (i in 3:length(phyloseq_list)) {
  glom_genus2 <- merge_phyloseq(glom_genus2, phyloseq_list[[i]])
}

# ---- Step 3: Perform pairwise ALDEx2 analysis (Wilcoxon tests) between location pairs ----
ListlOcaiton <- c('Control', 'EBPR', 'Vermont', 'AMF', 'EBPR_AMF', 'Vermont_AMF')
ListlOcaiton2 <- combn(ListlOcaiton, 2, simplify = FALSE)
wilcox.p_resDf <- NULL
for (cmb in ListlOcaiton2) {
  comparingdata <- subset_samples(glom_genus2, Location %in% cmb)
  print(cmb)
  
  # Define sample groupings for ALDEx2
  n1 <- sum(sample_data(comparingdata)$Location == cmb[1])
  n2 <- sum(sample_data(comparingdata)$Location == cmb[2])
  conds <- c(rep(cmb[1], n1), rep(cmb[2], n2))
  
  # Extract raw count matrix
  selex <- as.data.frame(otu_table(comparingdata))
  
  # Run ALDEx2 (Monte Carlo sampling = 128)
  x.all <- aldex(selex, conds, mc.samples = 128, test = "t",
                 effect = TRUE, include.sample.summary = FALSE,
                 denom = "all", verbose = FALSE, paired.test = FALSE)
  
  # ---- Step 4: Taxonomic table cleanup ----
  tax.clean <- as.data.frame(tax_table(comparingdata))
  tax.clean$Species <- rep("", nrow(tax.clean))
  tax.clean2 <- tax.clean
  for (i in 1:7) { tax.clean[, i] <- as.character(tax.clean[, i]) }
  tax.clean[is.na(tax.clean)] <- ""
  
  # Fill in missing taxonomy levels with placeholder
  for (i in 1:nrow(tax.clean)) {
    if (tax.clean[i, 2] == "") {
      kingdom <- paste("Kingdom:", tax.clean[i, 1], sep = "")
      tax.clean2[i, 2:7] <- kingdom
    } else if (tax.clean[i, 3] == "") {
      phylum <- paste("Phylum:", tax.clean[i, 2], sep = "")
      tax.clean2[i, 3:7] <- phylum
    } else if (tax.clean[i, 4] == "") {
      class <- paste("Class:", tax.clean[i, 3], sep = "")
      tax.clean2[i, 4:7] <- class
    } else if (tax.clean[i, 5] == "") {
      order <- paste("Order:", tax.clean[i, 4], sep = "")
      tax.clean2[i, 5:7] <- order
    } else if (tax.clean[i, 6] == "") {
      family <- paste("Family:", tax.clean[i, 5], sep = "")
      tax.clean2[i, 6:7] <- family
    } else if (tax.clean[i, 7] == "") {
      tax.clean2$Species[i] <- paste("Genus:", tax.clean$Genus[i], sep = "")
    }
  }
  
  # Build clean taxonomy reference
  List <- data.frame(
    rowname = rownames(tax.clean2),
    taxa = tax.clean2$Species,
    Phylum = tax.clean2$Phylum,
    Class = tax.clean2$Class,
    Order = tax.clean2$Order,
    Family = tax.clean2$Family,
    Genus = tax.clean2$Genus
  )
  
  # Merge ALDEx2 output with taxonomy
  x.all2 <- x.all
  x.all2$rowname <- rownames(x.all2)
  Output <- merge(x.all2, List, id = "rowname")
  A=gsub("rab.win.", "", colnames(Output)[3]);A
  B=gsub("rab.win.", "", colnames(Output)[4]);B

  # Extract significant results (BH-adjusted p-value < 0.05)

  res1<- Output %>% filter(wi.eBH < 0.05) %>%
    dplyr::select(rowname ,diff.btw, diff.win, effect, taxa, wi.eBH,
                  Phylum,Genus, Family,Order) %>%
    mutate(compare = ifelse(diff.btw > 0, 
                            paste(A, ">", B), 
                            paste(B, ">", A)))
  
  wilcox.p_resDf <- bind_rows(wilcox.p_resDf, res1)
}
# ---- Step 5: Filter strong-effect taxa for (effect size > 2.5) visualization ----

Data100 <- wilcox.p_resDf %>%
  filter(Absolute_effect >= 2.5) %>%
  subset(taxa != "Kingdom:Unassigned" & taxa != "Kingdom:Bacteria") %>%
  mutate(neg = -1 * diff.btw)

# ---- Step 6: Add specific taxa of interest (e.g. PAOs) ----
target_taxa <- c("Genus:Pseudomonas", "Genus:Thauera", "Genus:Paracoccus",
                 "Genus:Bacillus", "Genus:Rhodanobacter")
matched_taxa_df <- subset(wilcox.p_resDf, taxa %in% target_taxa) %>%
  subset(taxa != "Kingdom:Bacteria" & taxa != "Kingdom:Unassigned") %>%
  mutate(neg = -1 * diff.btw)

# Combine selected taxa for plotting
Data101 <- rbind(Data100, matched_taxa_df) %>% distinct()

# ---- Step 7: Normalize abundances and filter phyloseq object ----
otu_ids_to_keep <- unique(Data101$rowname)
glom_genus3 <- transform_sample_counts(glom_genus2, function(x) x / sum(x) * 100)
ref <- data.frame(otu_table(glom_genus2))

# Scale abundances across samples (centered and standardized)
transformed <- t(apply(data.frame(otu_table(glom_genus3)), 1, scale))
colnames(transformed) <- colnames(data.frame(otu_table(glom_genus3)))
rownames(transformed) <- rownames(data.frame(otu_table(glom_genus3)))
otu_table(glom_genus3) <- otu_table(as.matrix(transformed), taxa_are_rows = TRUE)

# ---- Step 8: Subset phyloseq object for selected taxa ----
glom_genus4 <- subset_taxa(glom_genus3, taxa_names(glom_genus3) %in% otu_ids_to_keep)

# ---- Step 9: Rewrite taxonomy with resolved species names ----
tax.clean <- as.data.frame(tax_table(glom_genus4))
tax.clean$Species <- rep("", nrow(tax.clean))
tax.clean2 <- tax.clean
for (i in 1:7) { tax.clean[, i] <- as.character(tax.clean[, i]) }
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)) {
  if (tax.clean[i, 2] == "") {
    kingdom <- paste("Kingdom:", tax.clean[i, 1], sep = "")
    tax.clean2[i, 2:7] <- kingdom
  } else if (tax.clean[i, 3] == "") {
    phylum <- paste("Phylum:", tax.clean[i, 2], sep = "")
    tax.clean2[i, 3:7] <- phylum
  } else if (tax.clean[i, 4] == "") {
    class <- paste("Class:", tax.clean[i, 3], sep = "")
    tax.clean2[i, 4:7] <- class
  } else if (tax.clean[i, 5] == "") {
    order <- paste("Order:", tax.clean[i, 4], sep = "")
    tax.clean2[i, 5:7] <- order
  } else if (tax.clean[i, 6] == "") {
    family <- paste("Family:", tax.clean[i, 5], sep = "")
    tax.clean2[i, 6:7] <- family
  } else if (tax.clean[i, 7] == "") {
    tax.clean2$Species[i] <- paste("Genus:", tax.clean$Genus[i], sep = "")
  }
}



# ---- Rewrite taxonomy table with cleaned hierarchy ----
tax.clean2 <- as.matrix(tax.clean2)  # Convert cleaned taxonomy to matrix format
rownames(tax.clean2) <- rownames(as.data.frame(tax_table(glom_genus4)))  # Preserve original rownames
head(tax.clean2)  # Preview first few entries for confirmation

# Update the phyloseq object with corrected taxonomy
tax_table(glom_genus4) <- tax_table(tax.clean2)

# ---- Prepare long-format abundance table for heatmap visualization ----
melted_genus4 <- psmelt(glom_genus4)  # Convert phyloseq object to dataframe (sample × taxon × abundance)

# ---- Determine number of unique genera for color palette sizing ----
num_colors <- length(unique(melted_genus4$Genus))
num_colors  # Print total to verify palette size

# ---- Create a color palette using Spectral gradient for distinct genus coloring ----
custom_palette <- colorRampPalette(brewer.pal(10, "Spectral"))(num_colors)

# ---- Reorder Location factor for consistent panel arrangement ----
melted_genus4$Location <- factor(melted_genus4$Location,
                                 levels = c("Control", "EBPR", "Vermont", "AMF", "EBPR_AMF", "Vermont_AMF"))

# ---- Arrange samples and remove blank genus labels ----
melted_genus4 <- melted_genus4 %>%
  arrange(Location) %>%                        # Sort by treatment group
  subset(Genus != "")                          # Drop entries without genus annotation

# ---- Set sample factor levels based on their order in the dataset ----
melted_genus4$Sample <- factor(melted_genus4$Sample, levels = unique(melted_genus4$Sample))

# ---- Remove entries where Species column is generic (e.g. starts with "Family:" or "Class:") ----
melted_genus4_filtered <- subset(melted_genus4, !grepl("^Family:|^Class:", Species))


dev.off()
# Create the heatmap
Aldex_Differential<-ggplot(melted_genus4_filtered, aes(x = Sample, y = Species, fill = Abundance)) +
  geom_tile(color = "white", width = 1, height = 1) +  # Adjusted height for tiles
  scale_fill_gradient2(
    low = "darkred", 
    mid = "white", 
    high = "darkblue", 
    midpoint = .02
  ) +
  theme_minimal() +
  theme_classic() +
  theme_bw() +  # Remove background
  labs(title = "Heatmap of Genus Abundance",
       x = "Sample",
       y = "Genus",
       fill = "Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.ticks.length = unit(0.2, "cm"),  # Adjust tick length
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +  # Adjust plot margins
  facet_grid(~ Location, scales = "free");Aldex_Differential  # Facet by location;


#Figure 4B
# Import genus-level relative abundance data formatted for microeco analysis
dataset <- microtable$new(otu_table = data.frame(otu_table(RA_data)), 
                          sample_table = data.frame(sample_data(RA_data)),
                          tax_table = data.frame(tax_table(RA_data)))
t1 <- trans_env$new(dataset = dataset, add_data = temp)
temp2<-temp %>% dplyr::select(Al, Ca, Fe, K, Mg, P, S, Zn)

t1 <- trans_env$new(dataset = dataset, add_data = temp2)
# Run Spearman correlation analyses
t1$cal_cor(
  use_data = "Species",             # Uses species-level microbial abundance data
  p_adjust_method = "BH",           # Applies Benjamini–Hochberg correction for multiple testing
  p_adjust_type = "Env",            # Adjusts p-values across environmental variables
  cor_method = "spearman"           # Computes non-parametric rank-based correlations
)
#Visualization
heatmap <- t1$plot_cor(
  filter_feature = c(""),
  pheatmap = TRUE,
  cluster_ggplot = "col",  
  color_palette = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))
); heatmap

#Figure 4C

# Import genus-level relative abundance data formatted for microeco analysis
temp2<-temp %>% dplyr::select(Al, Ca, Fe, K, Mg, P, S, Zn)
dataset2 <- microtable$new(otu_table = data.frame(otu_table(RA_data)), 
                           sample_table = data.frame(sample_data(RA_data)),
                           tax_table = data.frame(tax_table(RA_data)))
t1 <- trans_env$new(dataset = dataset2, add_data = temp2)
# Use add_abund_table parameter to add the extra data table (e.g. alpha diversity indices)
t1$cal_cor(add_abund_table = Vegan_diversity[,-4],p_adjust_method = "BH",
           p_adjust_type = "Env",cor_method='spearman')
# Visualization
Heatmapforalphadiversity<-t1$plot_cor();Heatmapforalphadiversity

# Figure 4D
# ---- Step 1: Convert KO abundance data into KEGG pathway abundance table ----
ko_abundance_file<-"KO.tsv"
abundance <- ko2kegg_abundance(ko_abundance_file)

# ---- Step 2: Define experimental groups ----
group <- "Group"
ListlOcaiton <- c('Control', 'EBPR', 'Vermont', 'AMF', 'AMF_EBPR', 'AMF_Vermont')

# ---- Generate all pairwise treatment combinations ----
ListlOcaiton2 <- combn(ListlOcaiton, 2, simplify = FALSE)

# ---- Prepare sample metadata ----
metadata <- data.frame(metadata)
metadata$Group <- factor(metadata$Group, 
                         levels = c("Control","EBPR","Vermont","AMF","AMF_EBPR","AMF_Vermont"))
rownames(metadata) <- metadata$sample_name

# ---- Step 3: Perform edgeR-based differential abundance analysis across treatment groups ----
# Note: This process is computationally intensive—proceed using the pre-saved 'Picrust_result' object for efficiency
Picrust_result <- NULL
for (treatment in ListlOcaiton2) {
  print(treatment)
  
  # Filter metadata to include samples for current comparison
  metadata2 <- metadata %>% filter(Group %in% treatment)
  metadata2$Group <- factor(metadata2$Group, levels = ListlOcaiton)
  
  # Subset abundance matrix for selected samples
  sample_name_variable <- metadata2$sample_name
  abundance2 <- abundance %>% dplyr::select(all_of(sample_name_variable))
  
  # Run pathway-level differential analysis via edgeR
  daa_results_df <- pathway_daa(abundance = abundance2,
                                metadata = metadata2,
                                group = group,
                                daa_method = "edgeR",
                                select = NULL,
                                p.adjust = "BH",
                                reference = NULL)
  
  # Annotate results with KEGG pathway info
  daa_results_df <- pathway_annotation(pathway = "KO",
                                       daa_results_df = daa_results_df,
                                       ko_to_kegg = TRUE)
  
  # ---- Step 4: Core edgeR modeling ----
  abundance_mat <- as.matrix(abundance2)
  dge <- edgeR::DGEList(counts = round(abundance_mat), group = metadata2$Group)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateCommonDisp(dge, verbose = TRUE)
  et <- edgeR::exactTest(dge, pair = c(1, 2))
  
  # Compile results into dataframe
  results <- data.frame(
    feature = rownames(abundance_mat),
    method = "edgeR",
    group1 = treatment[1],
    group2 = treatment[2],
    p_values = et$table$PValue,
    logFC = et$table$logFC,
    logCPM = et$table$logCPM,
    Comparison = rep(paste(treatment[1], 'vs', treatment[2]), nrow(et$table))
  )
  
  # Merge annotation and test statistics
  Merged <- merge(results, daa_results_df, by = "feature")
  Picrust_result <- bind_rows(Picrust_result, Merged)
}

# ---- Optionally reset to original result object (for quicker reuse) ----
Picrust_result <- Picrust_result_original

# ---- Step 5: Parse KEGG pathway hierarchy ----
Picrust_result$pathway_class2 <- sapply(strsplit(as.character(Picrust_result$pathway_class), ";"), `[`, 1)
Picrust_result$pathway_class3 <- sapply(strsplit(as.character(Picrust_result$pathway_class), ";"), `[`, 2)

# ---- Step 6: Filter significant functional pathways (non-disease-related) ----
summary_picrust_selected <- Picrust_result %>%
  filter(p_adjust < 0.05) %>%
  filter(pathway_class2 != "Human Diseases") %>%
  dplyr::select(feature, p_adjust, pathway_class, pathway_map,
                pathway_class2, pathway_class3, logFC, Comparison)

# ---- Step 7: Extract unique significant pathway features ----
SignificantFeature <- unique(summary_picrust_selected$feature)
length(SignificantFeature)
head(summary_picrust_selected)

# ---- Step 8: Build reference annotation for visualization ----
Reference <- summary_picrust_selected %>%
  dplyr::select(feature, pathway_map, pathway_class2, pathway_class3) %>% 
  distinct()

# ---- Step 9: Convert raw counts to relative abundance (%) ----
abundance100 <- abundance
abundance100_percentage <- sweep(abundance100, 2, colSums(abundance100), FUN = "/") * 100
abundance100_percentage$rownames <- rownames(abundance100_percentage)

# Subset to significant pathways only
Abundance_subset <- abundance100_percentage %>%
  filter(rownames %in% SignificantFeature)
Abundance_subset <- as.data.frame(Abundance_subset[,-ncol(Abundance_subset)])
Abundance_subset$feature <- rownames(Abundance_subset)

# Merge with annotation reference
Abundance_subset2 <- merge(Abundance_subset, Reference, by = "feature")

# ---- Step 10: Transpose abundance table to samples-by-feature format ----
abundance_t <- as.data.frame(t(Abundance_subset))
abundance_t$sample_name <- rownames(abundance_t)

# Merge with metadata
combined_data <- left_join(metadata, abundance_t, by = "sample_name")
rownames(combined_data) <- combined_data$sample_name

# ---- Step 11: Prepare abundance matrix for group-level summary ----
combined_data_fixed <- as.data.frame(combined_data[, -c(1, 3)])  # Drop unused metadata columns

# Convert KO columns to numeric
ko_numeric <- combined_data_fixed %>% 
  mutate(across(starts_with("ko"), ~ as.numeric(.)));ko_numeric

# Calculate group-wise average abundance
average_by_group <- ko_numeric %>%
  group_by(Group) %>%
  dplyr::summarise(across(starts_with("ko"), ~ mean(as.numeric(.x), na.rm = TRUE)))
average_by_group

# ---- Step 12: Filter features with ≥0.5% average abundance across sites ----
ko_sums <- colSums(average_by_group[, -1])
keep_cols <- names(ko_sums[ko_sums >= 0.5])
filtered_avg <- average_by_group[, c("Group", keep_cols)]
filtered_avg

# ---- Step 13: Z-score transformation ----
z_score_transformation <- function(df) {
  z_df <- as.data.frame(scale(df))
  return(z_df)
}
abundance100_z <- as.data.frame(z_score_transformation(
  as.data.frame(filtered_avg)[,-1]))  # Remove 'Group' before scaling

# ---- Step 14: Build composite data table with annotations ----
abundance100_z <- cbind(Group = average_by_group$Group, abundance100_z)
rownames(abundance100_z) <- abundance100_z$Group
abundance100_z_data <- abundance100_z[,-1]
COLNAMES <- colnames(abundance100_z_data)

# Row annotation (group identity)
row_annotation <- data.frame(Group = abundance100_z$Group)
rownames(row_annotation) <- rownames(abundance100_z)

# Column annotation (pathway class info)
col_annotation <- data.frame(
  pathway_class2 = summary_picrust_selected$pathway_class2[match(keep_cols, summary_picrust_selected$feature)],
  pathway_class3 = summary_picrust_selected$pathway_class3[match(keep_cols, summary_picrust_selected$feature)],
  pathway_map = summary_picrust_selected$pathway_map[match(keep_cols, summary_picrust_selected$feature)]
)
rownames(col_annotation) <- col_annotation$pathway_map
col_annotation <- col_annotation[, -3]  # Drop redundant column

# Align feature names with pathway maps
colnames(abundance100_z_data) <- summary_picrust_selected$pathway_map[match(keep_cols, summary_picrust_selected$feature)]
colnames(abundance100_z_data) == rownames(col_annotation)  # Ensure proper matching

# ---- Step 15: Setup color function and plot heatmap ----
col_fun <- colorRamp2(c(min(abundance100_z_data), 0, max(abundance100_z_data)), viridis(3))  # Blue-white-red

dev.off()  # Clear previous graphics device (if needed)

ht <- Heatmap(as.matrix(abundance100_z_data), 
              name = "Z-score",
              row_title = "Treatment",
              column_title = "Picrust2 defined function",
              row_names_gp = gpar(fontsize = 13, fontfamily = "Arial"),
              column_names_gp = gpar(fontsize = 13, fontfamily = "Arial"),
              top_annotation = HeatmapAnnotation(df = col_annotation),
              left_annotation = rowAnnotation(df = row_annotation),
              col = col_fun,
              column_order = order(col_annotation$pathway_class2),
              height = unit(3.2, "cm"),
              width = unit(20, "cm"),
              column_names_rot = 45); ht  # Render heatmap

save.image("C:/Users/noel1/OneDrive/Desktop/Cornell/Research/SARE/16analysis3/Natureplant/Rcodes/Fig4/Figure4_submission.RData")
