# Load necessary libraries
library(dplyr)
library(tidyr)
library(broom)
library(purrr)
library(multcompView)
library(biostat)
library(RColorBrewer)
library(ggplot2)
library(metagMisc)
library(phyloseq)
library(ggpubr)
library(reshape)
library(vegan)
library(plyr)    

################################################################
# ---- Load and Label Datasets for Analysis ----

# Load pre-processed phyloseq object containing ASV data and metadata
load("Fig3.RData")

# 'physeq3' holds the original (non-rarefied) ASV table, including taxa abundances and sample metadata
physeq3

# 'rarefied_phyloseq_original' contains the ASV dataset after rarefaction to standardize sequencing depth across samples
rarefied_phyloseq_original

#########################

# ---- Rarefaction: Equalizing Sampling Depth Across Samples ----

# Extract raw OTU table from physeq3 and convert to data frame
Rare_OTU <- data.frame(otu_table(physeq3))

# Transpose the table: rows → samples, columns → taxa
Rare_OTU_entire_t <- t(Rare_OTU)

# Compute observed species richness for each sample (pre-rarefaction)
S <- specnumber(Rare_OTU_entire_t)     # Richness vector
S                                       # Print for quick inspection

# Identify total reads per sample and minimum sample depth
raremax <- rowSums(Rare_OTU_entire_t)  # Total counts per sample
raremax                                # Distribution across samples
rarefying_number <- min(raremax)       # Smallest sample count

# Perform rarefaction to the smallest sample depth
rarefied_phyloseq <- rarefy_even_depth(physeq3,
                                       sample.size = rarefying_number,
                                       replace = FALSE,
                                       rngseed = 110)
rarefied_phyloseq<-rarefied_phyloseq_original

# Result: rarefied_phyloseq → normalized OTU abundance

# ---- Taxonomic Aggregation: Collapse to Genus-Level ----

# Aggregate OTUs to Genus level to reduce dimensionality and noise
Physeq3_ready <- rarefied_phyloseq %>%
  tax_glom(taxrank = "Genus", NArm = FALSE)

# ---- Prevalence Filtering by Location ----

# Initialize list to store filtered phyloseq objects per location
phyloseq_list <- list()

# Loop through each unique location in metadata
for (i in unique(sample_data(Physeq3_ready)$Location)) {
  print(i)  # Track progress
  
  # Subset samples for current location
  Individualdata <- subset_samples(Physeq3_ready, Location %in% i)
  
  # Filter genera: retain those present in ≥70% of samples
  tmp_data <- phyloseq_filter_prevalence(Individualdata, prev.trh = 0.70)
  
  # Save location-specific filtered phyloseq object
  phyloseq_list[[i]] <- tmp_data
}

# ---- Merge Location-Specific Filtered Phyloseq Objects ----

# Initialize empty phyloseq object
combined_site <- NULL

# Merge all location-wise filtered data (assumes 6 locations)
for (i in 1:6) {
  print(paste("Treatment is", i))
  print(phyloseq_list[[i]])
  
  # Merge current location data with cumulative object
  combined_site <- merge_phyloseq(combined_site, phyloseq_list[[i]])
}

# Final output: combined_site = ready for analysis (e.g. SPIEC-EASI, diversity, ordination)

# ---- Filter taxa with total counts >100 across all samples ----
Morethan100 <- prune_taxa(taxa_sums(Physeq3_ready) > 100, Physeq3_ready)

# ---- Transform raw counts to relative abundance (%) per sample ----
Morethan100_1 <- transform_sample_counts(Morethan100, function(x) 100 * x / sum(x))

# ---- Collapse OTU table to Phylum level ----
glom_phylum <- tax_glom(Morethan100_1, taxrank = 'Phylum', NArm = FALSE)

# ---- Convert phyloseq object into long format data frame for plotting ----
melted_phylum <- psmelt(glom_phylum)

# ---- Generate color palette for number of unique phyla ----
num_colors <- length(unique(melted_phylum$Phylum))
num_colors
custom_palette <- colorRampPalette(brewer.pal(10, "Spectral"))(num_colors)

# ---- Re-melt in case of upstream transformation changes ----
melted_phylum <- psmelt(glom_phylum)

# ---- Reorder Location factor for consistent faceting ----
melted_phylum$Location <- factor(melted_phylum$Location,
                                 levels = c("Control", "EBPR", "Vermont", "AMF", "EBPR_AMF", "Vermont_AMF"))

# ---- Arrange samples by Location block for plot clarity ----
melted_phylum2 <- melted_phylum %>% arrange(Location)
melted_phylum2$Sample <- factor(melted_phylum2$Sample, levels = unique(melted_phylum2$Sample))

# ---- Create stacked bar plot at Phylum level ----
Phylum_distribution <- ggplot(melted_phylum2, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) + # Bar width and stacking
  scale_fill_manual(values = custom_palette) + # Custom color palette
  theme_minimal() +
  theme_classic() +
  theme_bw() +     # Remove grey background
  labs(title = "Compacted Stacked Bar Graph",
       x = "Sample",
       y = "Abundance",
       fill = "Phylum") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  facet_grid(~ Location, scale = "free")     # Facet by Location
Phylum_distribution

# ---- Collapse OTU table to Genus level ----
glom_genus <- tax_glom(Morethan100_1, taxrank = 'Genus', NArm = FALSE)

# ---- Transform to relative abundance again for genus filtering ----
glom_genus2 <- transform_sample_counts(glom_genus, function(x) x / sum(x) * 100)

# ---- Filter genera with mean abundance > 0.5% ----
glom_genus3 <- filter_taxa(glom_genus2, function(x) mean(x) > 0.5, TRUE)

# ---- Melt phyloseq to long format for Genus-level plotting ----
melted_genus <- psmelt(glom_genus3)

# ---- Generate color palette based on number of unique genera ----
num_colors <- length(unique(melted_genus$Genus))
num_colors
custom_palette <- colorRampPalette(brewer.pal(10, "Spectral"))(num_colors)

# ---- Reorder Location and Sample factors for plot consistency ----
melted_genus$Location <- factor(melted_genus$Location,
                                levels = c("Control", "EBPR", "Vermont", "AMF", "EBPR_AMF", "Vermont_AMF"))
melted_genus2 <- melted_genus %>% arrange(Location)
melted_genus2$Sample <- factor(melted_genus2$Sample, levels = unique(melted_genus2$Sample))

# ---- Create stacked bar plot at Genus level ----
Genus_distribution <- ggplot(melted_genus2, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  scale_fill_manual(values = custom_palette) +
  theme_minimal() +
  theme_classic() +
  theme_bw() +
  labs(title = "Compacted Stacked Bar Graph",
       x = "Sample",
       y = "Abundance",
       fill = "Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  facet_grid(~ Location, scale = "free")
Genus_distribution

#Fig 3B
# ---- Aggregate OTU table to Phylum level ----
glom2 <- tax_glom(Morethan100_1, taxrank = 'Phylum', NArm = FALSE)

# ---- Convert phyloseq object to long format dataframe for plotting and stats ----
ps.melt3 <- psmelt(glom2)

# ---- Filter phyla based on relative abundance: retain phyla with median > 0.5% ----
medians <- ddply(ps.melt3, ~Phylum, function(x) c(median = median(x$Abundance)))  # Compute median per phylum
notselected <- medians[round(medians$median, digits = 1) <= 0.5, ]$Phylum         # Phyla to combine into "<0.5%" group
selected <- medians[round(medians$median, digits = 1) > 0.5, ]$Phylum              # Phyla to keep separately
selected  # Preview selected phyla

# ---- Rename low-abundance phyla for plot simplification ----
ps.melt3[ps.melt3$Phylum %in% notselected, ]$Phylum <- "Phyla < 0.5% abundance"
ps.melt3  # Updated dataframe

# ---- Run Kruskal-Wallis test for each selected phylum across locations ----
phyla_significance <- ps.melt3 %>%
  filter(Phylum != "Phyla < 0.5% abundance") %>%
  nest(data = -Phylum) %>%
  mutate(test = map(data, ~kruskal.test(Abundance ~ Location, data = .x) %>% tidy)) %>%
  unnest(test) %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH")) %>%
  dplyr::select(Phylum, p.adjust)

# ---- Perform pairwise Wilcoxon tests with BH adjustment ----
pairwise_comparisons <- ps.melt3 %>%
  filter(Phylum != "Phyla < 0.5% abundance") %>%
  group_by(Phylum) %>%
  nest() %>%
  mutate(pairwise_tests = map(data, ~pairwise.wilcox.test(.x$Abundance, .x$Location,
                                                          p.adjust.method = "BH") %>% tidy)) %>%
  unnest(pairwise_tests) %>%
  dplyr::select(Phylum, group1, group2, p.value) %>%
  group_by(Phylum) %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH")) %>%
  ungroup()
pairwise_comparisons  # Display results

# ---- Extract list of statistically relevant phyla based on pairwise tests ----
selected_phyla <- sort(unique(pairwise_comparisons$Phylum))

# ---- Generate compact letter display for each phylum using pairwise tests ----
get_letters <- function(df) {
  comparisons <- with(df, pairwise.wilcox.test(Abundance, Location, p.adjust.method = "BH"))
  letters <- multcompLetters(comparisons$p.value)$Letters
  return(letters)
}

significance_letters <- ps.melt3 %>%
  group_by(Phylum) %>%
  nest() %>%
  mutate(letters = map(data, get_letters)) %>%
  unnest(data) %>%
  unnest(letters)
head(significance_letters)  # Preview letter assignments

# ---- Subset phyloseq object to retain only selected phyla ----
glom2_subset <- subset_taxa(glom2, Phylum %in% selected_phyla)

# ---- Convert subset phyloseq to long format for plotting ----
glom2_melt <- psmelt(glom2_subset)

# ---- Initialize plot list and significance table containers ----
GGPLOT_LIST <- list()
significance_table <- NULL

# ---- Generate boxplots for each selected phylum ----
for (phylum in selected_phyla) {
  glom2_melt2 <- glom2_melt %>% filter(Phylum == phylum)  # Filter data for current phylum
  print(phylum)
  glom2_melt2_z <- glom2_melt2 %>%
    mutate(Abundance_z = scale(Abundance))
  
  data333<-pairwise.wilcox.test(glom2_melt2_z$Abundance ,  
                                glom2_melt2_z$Location);data333
  Significance_letters<- biostat::make_cld(data333);Significance_letters
  colnames(Significance_letters)[1]<-"Location";Significance_letters
  Significance_letters2<-Significance_letters
  Significance_letters2$Phylum<-rep(phylum,nrow(Significance_letters2))

    # ---- Build ggplot boxplot for phylum across locations ----
  Phylum_plot <- ggplot(glom2_melt2, aes(x = Location, y = Abundance, fill = Location)) +
    geom_boxplot(aes(color = Location), outlier.size = 0.4, width = 0.8, alpha = 0.5,
                 position = position_dodge(width = 0.9)) +
    geom_point(aes(color = Location), position = position_jitterdodge(dodge.width = 0.9)) +
    theme_bw() +
    geom_text(data = Significance_letters,aes( x = Location , 
                                               y = max(glom2_melt2$Abundance), label = cld) ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_manual(values = c('seagreen3', 'cornflowerblue', 'orange2', 'coral',
                                 'royalblue4', 'coral4')) +
    scale_color_manual(values = c('seagreen3', 'cornflowerblue', 'orange2', 'coral',
                                  'royalblue4', 'coral4')) +
    labs(title = phylum, y = "Relative abundance (%)") +
    theme_classic() +
    theme(legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 12),
          axis.title.x = element_blank(),
          legend.position = "none",
          axis.title.y = element_text(size = 12, face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          aspect.ratio = 0.8)
  Phylum_plot
  
  # ---- Save to plot list ----
  GGPLOT_LIST[[phylum]] <- Phylum_plot
}

# ---- Define layout matrix for multi-panel display ----
lay <- rbind(c(1, 2),
             c(3, 4),
             c(5, 6),
             c(7, 8))

# ---- Arrange phylum boxplots in a 4x2 grid ----
grid.arrange(GGPLOT_LIST[[1]], GGPLOT_LIST[[2]], GGPLOT_LIST[[3]], GGPLOT_LIST[[4]],
             GGPLOT_LIST[[5]], GGPLOT_LIST[[6]], GGPLOT_LIST[[7]], GGPLOT_LIST[[8]],
             layout_matrix = lay, respect = TRUE)


#Fig 3C
# ---- Calculate relative abundance (%), excluding ambiguous genus annotations ----
RA_data <- combined_site %>%
  transform_sample_counts(function(x) {x / sum(x) * 100}) %>%
  subset_taxa(Genus != "-1")  # Remove taxa with unidentified genus

# ---- Run overall PERMANOVA (Adonis2) to test treatment effect on microbial composition ----
Matrix <- t(as.matrix(otu_table(RA_data)))                      # Transpose OTU table
metadata <- as(sample_data(RA_data), "data.frame")              # Extract metadata
adonis_data <- adonis2(Matrix ~ Location, data = metadata)      # Run PERMANOVA
adonis_data                                                     # View results

# ---- Compute Bray–Curtis dissimilarities and perform NMDS ordination ----
iDist <- distance(RA_data, method = "bray")                      # Bray–Curtis distance
iMDS <- ordinate(RA_data, "NMDS", distance = iDist)              # Non-metric multidimensional scaling
iMDS                                                             # NMDS ordination object

# ---- Fit continuous environmental variables (e.g. mineral concentrations) to NMDS ordination ----
en <- envfit(iMDS, temp2, permutations = 999, na.rm = TRUE)      # Spearman correlation test for vectors
en                                                               # Print result table
'''
***VECTORS

      NMDS1    NMDS2     r2 Pr(>r)   
Al -0.88799 -0.45986 0.0033  0.910   
Ca  0.78518  0.61927 0.0135  0.657   
Fe -0.31165  0.95020 0.0205  0.551   
K   0.21635 -0.97632 0.2045  0.002 **
Mg  0.66064 -0.75070 0.1014  0.045 * 
P   0.98665  0.16284 0.1446  0.013 * 
S   0.97681  0.21410 0.1260  0.028 * 
Zn  0.79969 -0.60041 0.0098  0.757   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999
'''
# ---- Prepare arrow coordinates for environmental vectors ----
en_coord_cont <- data.frame(
  NMDS1 = en$vectors$arrows[1:8],
  NMDS2 = en$vectors$arrows[9:16]
)
rownames(en_coord_cont) <- c('Al','Ca','Fe','K','Mg','P','S','Zn')

# ---- Extract site scores from ordination ----
NMDS <- data.frame(
  MDS1 = iMDS$points[,1],
  MDS2 = iMDS$points[,2],
  group = sample_data(RA_data)$Location
)

# ---- Compute mean coordinates by group for ellipse centers ----
NMDS.mean <- aggregate(NMDS[,1:2], list(group = NMDS$group), mean)

# ---- Define ellipse generation function (confidence boundary based on covariance) ----
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# ---- Construct ellipse coordinates for each location group ----
df_ell <- data.frame()
for (g in levels(NMDS$group)) {
  df_ell <- rbind(df_ell, cbind(
    as.data.frame(with(
      NMDS[NMDS$group == g,],
      veganCovEllipse(cov.wt(cbind(MDS1, MDS2), wt = rep(1/length(MDS1), length(MDS1)))$cov,
                      center = c(mean(MDS1), mean(MDS2)))
    )),
    group = g
  ))
}
df_ell$group <- factor(df_ell$group,
                       levels = c("Control", "EBPR", "Vermont", "AMF", "EBPR_AMF", "Vermont_AMF"))

# ---- Rescale environmental vector lengths for visualization ----
en_coord_cont2 <- 0.6 * en_coord_cont

# ---- Build NMDS scatter plot with ellipses and environmental vectors ----
NMDS_Sorghum2 <- ggplot(data = NMDS, aes(MDS1, MDS2)) +
  geom_point(aes(color = group), size = 0.8) +
  geom_path(data = df_ell, aes(x = MDS1, y = MDS2, colour = group),
            size = 0.6, linetype = 2, alpha = 1) +
  scale_shape_manual(values = c(16,17,18,15,25,8)) +
  scale_color_manual(values = c('seagreen3','cornflowerblue','orange2',
                                'coral','royalblue4','coral4')) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_segment(data = en_coord_cont2,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               size = 1, alpha = 0.5, colour = "cornsilk3",
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = en_coord_cont2, aes(x = NMDS1, y = NMDS2),
            colour = "cornsilk4", fontface = "bold",
            label = row.names(en_coord_cont2))
NMDS_Sorghum2  # Display plot

# ---- Define color scheme for Location groups ----
Location_c <- c('seagreen3','cornflowerblue','orange2','coral','royalblue4','coral4')

# ---- MDS1 boxplot per treatment group ----
xplot <- ggplot(data = NMDS, aes(x = group, y = MDS1, color = group, fill = group)) +
  geom_boxplot(alpha = 0.5, size = 0.5, width = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.5), shape = 21) +
  theme_classic() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = Location_c) +
  scale_color_manual(values = Location_c) +
  clean_theme() + rotate() +
  theme(aspect.ratio = 0.12, legend.position = "none")
xplot  # Display plot

# ---- MDS2 boxplot per treatment group ----
yplot <- ggplot(data = NMDS, aes(x = group, y = MDS2, color = group, fill = group)) +
  geom_boxplot(alpha = 0.5, size = 0.5, width = 0.9, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.5), shape = 21) +
  theme_classic() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = Location_c) +
  scale_color_manual(values = Location_c) +
  clean_theme() +
  theme(aspect.ratio = 6, legend.position = "none")
yplot  # Display plot

# ---- Loop over all pairwise group comparisons for PERMANOVA ----
wilcox.p_resDf <- NULL
for (cmb in list(
  c("Control", "EBPR"), c("Control", "Vermont"), c("Control", "AMF"),
  c("Control", "EBPR_AMF"), c("Control", "Vermont_AMF"), c("EBPR", "Vermont"),
  c("EBPR", "AMF"), c("EBPR", "EBPR_AMF"), c("EBPR", "Vermont_AMF"),
  c("Vermont", "AMF"), c("Vermont", "EBPR_AMF"), c("Vermont", "Vermont_AMF"),
  c("AMF", "EBPR_AMF"), c("AMF", "Vermont_AMF"), c("EBPR_AMF", "Vermont_AMF")
)) {
  print(cmb)
  comparingdata <- phyloseq::subset_samples(RA_data, Location %in% cmb)
  Matrix <- t(as.matrix(otu_table(comparingdata)))
  metadata <- as(sample_data(comparingdata), "data.frame")
  adonis_data <- adonis2(Matrix ~ Location, data = metadata)
  
  res <- adonis_data %>%
    summarise(Rvalue = adonis_data$R2[1]) %>%
    mutate(Pvalue = adonis_data$'Pr(>F)'[1]) %>%
    mutate(Comparison = paste(cmb, collapse = " vs "))
  wilcox.p_resDf <- bind_rows(wilcox.p_resDf, res)
}

# ---- Apply Benjamini–Hochberg correction for multiple comparisons ----
wilcox.p_resDf$BHmodifiedP <- p.adjust(wilcox.p_resDf$Pvalue, method = "BH")
wilcox.p_resDf$number <- 1:15
# ---- Assign custom hex colors for each pairwise comparison ----
wilcox.p_resDf$color<-c("3E7DAFFF", "458BB3FF","6BD7C9FF","4081B0FF","7EFCD3FF",
                        "04078DFF","3A73ACFF","5EBCC1FF","4E9DB8FF","2850A2FF",
                        "00008BFF","5AB6BFFF","4285B1FF","6CD9C9FF","4E9DB8FF")
Rdata<-wilcox.p_resDf

# ---- Create bar plot of Adonis R² values across pairwise group comparisons ----

Rplot<- ggplot(data=Rdata,aes(x=reorder(Comparison, number,decreasing=TRUE),y=Rvalue)) +
  geom_bar(aes(fill =Rvalue ),alpha=0.8,width=0.9,stat="identity",show.legend = FALSE) +
  theme_classic() +
  theme_bw()+ #Remove background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_gradient(low="cadetblue2", high="blue4") +
  theme(axis.text.x = element_text( size=9, angle=0))+
  theme(axis.text.y = element_text( size=9, angle=0))+
  xlab("Adonis R2")+
  theme(aspect.ratio = 4)+
  scale_y_continuous(expand = c(0, 0.02, 0, 0.25))+
  coord_flip();Rplot






#Fig 3D
# ---- Convert OTU table to sample-by-taxon matrix ----
datasetforhill <- t(as.matrix(otu_table(RA_data)))  # Transpose so samples are rows
head(datasetforhill)                                # Preview first few rows

# ---- Calculate alpha diversity metrics ----
H <- diversity(datasetforhill)           # Shannon diversity index (entropy-based)
S <- specnumber(datasetforhill)          # Species richness (number of taxa per sample)
J <- H / log(S); J                       # Pielou’s evenness: H normalized by log(species richness)

# ---- Store diversity metrics in one data frame ----
Vegan_diversity <- data.frame(H, S, J)               # Columns: Shannon, richness, evenness
Vegan_diversity$Location <- sample_data(physeq3)$Location  # Append sample metadata

# ---- Reshape the data for ggplot: long format with 'variable' and 'value' columns ----
Vegan_diversity2 <- melt(Vegan_diversity, id = c('Location'))

# ---- Build boxplots for each diversity metric across Locations ----
withhill <- ggplot(data = Vegan_diversity2, aes(x = Location, y = value,
                                                fill = Location, color = Location)) +
  facet_wrap(~variable, scales = "free") +                  # One subplot per metric
  
  # ---- Boxplot visualization with jittered points overlaid ----
geom_boxplot(outlier.size = 0.4, width = 0.8, alpha = 0.5,
             position = position_dodge(width = 0.9)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9)) +
  
  # ---- Apply visual themes to clean the background ----
theme_minimal() +
  theme_classic() +
  theme_bw() +                                              # Force white background
  
  # ---- Remove grid lines for clarity ----
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
  
  # ---- Apply custom color palettes to treatments ----
scale_color_manual(values = c('seagreen3', 'cornflowerblue', 'orange2',
                              'coral', 'royalblue4', 'coral4')) +
  scale_fill_manual(values = c('seagreen3', 'cornflowerblue', 'orange2',
                               'coral', 'royalblue4', 'coral4')) +
  
  # ---- Customize plot borders and axis appearance ----
theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +                   # Hide x-axis text and ticks
  
  theme(axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12, face = "bold")) +
  
  # ---- Set fixed aspect ratio and y-axis expansion ----
theme(aspect.ratio = 1) +
  scale_y_continuous(expand = c(0.1, 0, 0.15, 0))           # Pad axis boundaries

# ---- Render the final faceted plot of alpha diversity metrics ----
withhill
