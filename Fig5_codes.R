
# ---- Load Libraries for Fig.5 Network Analysis ----
library(SpiecEasi)
library(phyloseq)
library(Matrix)
library(pulsar)
library(batchtools)
library(ggrepel)
library(RCy3)
library(igraph)
library(dplyr)
library(plyr)
library(htmlwidgets)
library(webshot2)
library(tidygraph)
library(RColorBrewer)
library(chorddiag)


# ---- Loaded Analysis Datasets ----
load("Fig5.RData")
# --- Run SPIEC-EASI network construction --- 
# This loop executes a computationally intensive SPIEC-EASI pipeline.
# By allocating 100 cores and setting rep.num = 1000, the stability of 
# microbial association network estimates is significantly enhanced.
# However, such parallelization places a considerable load on high-performance systems,
# and may require runtime monitoring or checkpointing in cluster environments.

# The output of this process is stored as 'Network_saved2_original',
# which contains the previously generated SPIEC-EASI network objects.

pargs2 <- list(rep.num=1000, seed=10010, ncores=100)

# Loop over each unique location in the sample metadata
for (i in unique(sample_data(combined_site)$Location)) {
  print(i)  # Print the current location name for monitoring progress
  # Subset the phyloseq object to include only samples from the current location
  ps <- phyloseq::subset_samples(combined_site, Location %in% i)
  # Prune taxa with zero total counts (i.e., remove taxa absent in all samples for this location)
  amgut2.filt.phy <- prune_taxa(taxa_sums(ps) > 0, ps)
  # Estimate microbial association network using the Meinshausen-Bühlmann method 
  # SPIEC-EASI parameters:
  # - lambda.min.ratio: Minimum lambda threshold for regularization path
  # - nlambda: Number of lambda values to evaluate
  # - sel.criterion: Model selection criterion; here, 'stars' uses stability-based selection
  # - pulsar.select: Whether to perform model selection with Pulsar backend
  # - pulsar.params: Additional Pulsar tuning parameters provided via `pargs2`
  se.mb.amgut2 <- spiec.easi(amgut2.filt.phy,
                             method = 'mb',
                             lambda.min.ratio = 1e-4,
                             nlambda = 100,
                             sel.criterion = 'stars',
                             pulsar.select = TRUE,
                             pulsar.params = pargs2)
  # Convert the adjacency matrix from SPIEC-EASI output into an igraph object
  # - `getRefit()` extracts the optimal network
  # - `vertex.attr` assigns taxa names as node labels
  ig2.mb <- adj2igraph(getRefit(se.mb.amgut2),
                       vertex.attr = list(name = taxa_names(amgut2.filt.phy)))
  # Save the network estimation results, igraph object, and filtered phyloseq object
  # into a named list indexed by the location
  Network_saved2[[i]] <- list(
    se.mb.amgut2 = se.mb.amgut2,     # SPIEC-EASI model object
    ig2.mb = ig2.mb,                 # Microbial interaction network as igraph
    amgut2.filt.phy = amgut2.filt.phy  # Preprocessed phyloseq object for this location
  )
}

# ---- Restore Recommended Network Object ----
# Use the original version of 'Network_saved2' for downstream analyses
Network_saved2 <- Network_saved2_original

list_hub<- list()

for (treatment in names(Network_saved2)) {
  print(treatment)  # Track progress by printing treatment name
  
  # Extract SPIEC-EASI model and filtered phyloseq object
  se.mb.amgut2 <- Network_saved2[[treatment]][["se.mb.amgut2"]]
  amgut2.filt.phy <- Network_saved2[[treatment]][["amgut2.filt.phy"]]
  
  # Symmetrize inverse covariance matrix (beta) to get edge weights
  optbeta <- as.matrix(symBeta(getOptBeta(se.mb.amgut2)))
  
  # Assign edge colors: blue for positive, red for negative associations
  edge_cols <- ifelse(optbeta > 0, 'darkblue', 'darkred')[upper.tri(optbeta) & optbeta != 0]
  
  # Convert adjacency matrix to igraph object, with taxa names and edge colors
  ig2.mb <- adj2igraph(getRefit(se.mb.amgut2), rmEmptyNodes = TRUE,
                       vertex.attr = list(name = taxa_names(amgut2.filt.phy)),
                       edge.attr = list(color = edge_cols))
  
  # Add node degree (number of connections) to igraph vertices
  V(ig2.mb)$degree <- igraph::degree(ig2.mb)
  
  # Retrieve OTU IDs used in regression-based covariance estimation
  otu.ids <- colnames(se.mb.amgut2[[1]]$data)
  
  # Manually re-derive edge colors and weights for finer control
  betaMat <- as.matrix(symBeta(getOptBeta(se.mb.amgut2)))
  edge.colors <- c()
  edge.weight <- c()
  
  for (e.index in seq_along(E(ig2.mb))) {
    adj.nodes <- ends(ig2.mb, E(ig2.mb)[e.index])
    xindex <- which(otu.ids == adj.nodes[1])
    yindex <- which(otu.ids == adj.nodes[2])
    beta <- betaMat[xindex, yindex]
    
    weight <- format(beta, digits = 3)
    edge.weight <- append(edge.weight, weight)
    edge.colors <- append(edge.colors, ifelse(beta > 0, "darkblue", "red"))
  }
  
  # Apply computed attributes to igraph
  E(ig2.mb)$color <- edge.colors
  E(ig2.mb)$weight <- as.numeric(edge.weight)
  E(ig2.mb)$weight_abs <- abs(E(ig2.mb)$weight)
  V(ig2.mb)$number <- seq_along(V(ig2.mb))
  
  # Remove zero-degree vertices and extract positive/negative subgraphs
  ig2.mb2 <- delete_vertices(ig2.mb, V(ig2.mb)$degree == 0)
  g_70_positive <- delete_edges(ig2.mb2, E(ig2.mb2)[E(ig2.mb2)$weight < 0])
  g_70_negative <- delete_edges(ig2.mb2, E(ig2.mb2)[E(ig2.mb2)$weight > 0])
  
  # Convert igraph to data frames for node and edge attributes
  vertices_70 <- igraph::as_data_frame(ig2.mb2, what = "vertices")
  relations_70 <- igraph::as_data_frame(ig2.mb2, what = "edges")
  
  # Reconstruct igraph from raw data frames (more flexible for manipulation)
  igraph_g70 <- graph_from_data_frame(relations_70, directed = FALSE, vertices = vertices_70)
  
  # Compute betweenness centrality and detect modules via edge betweenness clustering
  betweennessCentrality <- betweenness(g_70_positive)
  mst.communities <- edge.betweenness.community(g_70_positive)
  mst.clustering <- make_clusters(g_70_positive, membership = mst.communities$membership)
  
  # Annotate graph nodes with community, degree, and centrality metrics
  V(g_70_positive)$membership <- mst.communities$membership
  V(g_70_positive)$degree <- degree(g_70_positive)
  V(g_70_positive)$BetweennessCentrality <- betweennessCentrality
  E(g_70_positive)$edge.betweenness <- mst.communities$edge.betweenness
  
  # Calculate Fruchterman-Reingold layout coordinates
  coords <- layout_with_fr(g_70_positive)
  V(g_70_positive)$layout_first <- coords[, 1]
  V(g_70_positive)$layout_second <- coords[, 2]
  
  # Prepare graph for ggplot by joining node metadata
  Org_node <- data.frame(
    ASV = V(g_70_positive)$name,
    BetweennessCentrality = V(g_70_positive)$BetweennessCentrality,
    degree = V(g_70_positive)$degree,
    membership = V(g_70_positive)$membership
  )
  
  # Normalize network centrality metrics and identify top-performing hub taxa in the upper 10% 
  #by both betweenness centrality and node degree
  standarizeddata <- sqrt(Org_node[, 2:3])
  BC10 <- 0.9 * max(standarizeddata$BetweennessCentrality)
  Degree10 <- 0.9 * max(standarizeddata$degree)
  
  combineddata <- cbind(Org_node[, 1], standarizeddata, Org_node[, 4])
  colnames(combineddata) <- colnames(Org_node)
  
  # Label hubs if they pass thresholds for both degree and centrality
  df <- combineddata[, 1:3]
  df$color <- ifelse(df$BetweennessCentrality >= quantile(df$BetweennessCentrality, 0.9) &
                       df$degree >= quantile(df$degree, 0.9), "darkred", "darkgray")
  a_select <- df[df$color == "darkred", ]$ASV
  
  # Merge taxonomy to hub table for visualization
  taxtable100 <- data.frame(tax_table(amgut2.filt.phy))
  taxtable100$ASV <- rownames(taxtable100)
  df2 <- cbind(Org_node, df[, "color"])
  colnames(df2)[5] <- "color"
  ggplottable <- merge(df2, taxtable100, by = 'ASV')
  
  # Generate ggplot with hub species annotated
  ggplot_hub <- ggplot(ggplottable, aes(x = degree, y = BetweennessCentrality)) +
    geom_point(size = 3, aes(color = color)) +
    geom_text_repel(data = subset(ggplottable, color == "darkred"), aes(label = Species), show.legend = FALSE) +
    scale_color_manual(values = c("darkgrey", "darkred")) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_text(size = 14),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14))
  
  # Convert the positive subnetwork into vertex and edge data frames
  
  vertices_positive<-as_data_frame(g_70_positive, what="vertices")
  relations_positive<-as_data_frame(g_70_positive, what="edges")
  
  # Rebuild the positive association igraph object using data frames (undirected network)
  
  igraph_positive <- graph_from_data_frame(relations_positive, directed=FALSE, 
                                           vertices=vertices_positive)
  
  # Extract subnetwork around hub taxa
  Conv_pos <- igraph_positive
  Conv_workingtable2 <- ggplottable %>% filter(color == "darkred") %>% arrange(color)
  Conv_hubspecies <- Conv_workingtable2$ASV
  Conv_hub_network <- neighbors(Conv_pos, Conv_hubspecies, mode = "total")
  Conv_hub_network_igraph <- induced_subgraph(Conv_pos, c(Conv_hub_network, V(Conv_pos)[Conv_hubspecies]))
  
  # Save all processed objects to output list
  list_hub[[treatment]] <- list(
    OriginalIgraph = Conv_pos,
    a_select = a_select,
    ig2.mb = ig2.mb,
    ig2.mb2 = ig2.mb2,
    g_70_positive = g_70_positive,
    g_70_negative = g_70_negative,
    vertices_positive = vertices_positive,
    relations_positive = relations_positive,
    igraph_positive = igraph_positive,
    ggplottable = ggplottable,
    ggplot_hub = ggplot_hub,
    Hubnetwork_igaph = Conv_hub_network_igraph
  )
}

# Create combined hub networks

# Initialize storage objects
Hub_relations_g_union <- NULL
Hub_phyloseq_union <- NULL
ListlOcaiton<-c('Control', 'EBPR', 'Vermont' ,'AMF' ,'EBPR_AMF', 'Vermont_AMF');ListlOcaiton

# Loop through each location/treatment
for (a in ListlOcaiton) {
  print(a)  # Monitor progress
  
  # Extract hub subnetwork for current treatment
  Conv_hub_network_igraph <- list_hub[[a]][["Hubnetwork_igaph"]]
  
  # Convert igraph to edge dataframe and tag treatment
  added <- igraph::as_data_frame(Conv_hub_network_igraph, what = "edges")
  added$Treatment <- rep(a, nrow(added))
  
  # Combine edge data across treatments
  Hub_relations_g_union <- rbind(Hub_relations_g_union, added)
  
  # Merge phyloseq objects across treatments for annotation
  amgut2.filt.phy <- Network_saved2[[a]][["amgut2.filt.phy"]]
  Hub_phyloseq_union <- merge_phyloseq(Hub_phyloseq_union, amgut2.filt.phy)
}
# Retrieve taxonomy
taxtable100 <- data.frame(tax_table(Hub_phyloseq_union))
# Reformat for source node ('from') annotation
from_taxtable100 <- taxtable100
from_taxtable100$from <- rownames(from_taxtable100)
colnames(from_taxtable100)[1:7] <- c("from_Kingdom", "from_Phylum", "from_Class",
                                     "from_Order", "from_Family", "from_Genus", "from_Species")
# Reformat for target node ('to') annotation
to_taxtable100 <- taxtable100
to_taxtable100$to <- rownames(to_taxtable100)
colnames(to_taxtable100)[1:7] <- c("to_Kingdom", "to_Phylum", "to_Class",
                                   "to_Order", "to_Family", "to_Genus", "to_Species")
# Merge taxonomy with edge data
Hub_relations_g_union2 <- merge(Hub_relations_g_union, from_taxtable100, by = "from")
Hub_relations_g_union3 <- merge(Hub_relations_g_union2, to_taxtable100, by = "to") 

# Hub_relations_g_union3 : Ready for Cytoscape import and comparative visualization!



#Fig. 5B

# Loop through each treatment
for (treatment in names(list_hub)) {
  print(treatment)
  # Extract full igraph and convert to edge list
  g_positive <- list_hub[[treatment]][["ig2.mb"]]
  vertices_pos <- igraph::as_data_frame(g_positive, what = "vertices")
  relations_pos <- igraph::as_data_frame(g_positive, what = "edges")
  # Set working directory and load taxonomic metadata
  taxa22 <- data.frame(tax_table(combined_site))
  taxa22$ASV <- rownames(taxa22)
  # Reformat edge endpoints
  source <- relations_pos["from"]; colnames(source) <- "ASV"
  target <- relations_pos["to"]; colnames(target) <- "ASV"
  # Merge taxonomy to source and target OTUs
  Output_source <- merge(source, taxa22, by = "ASV")
  Output_target <- merge(target, taxa22, by = "ASV")
  # Extract Phylum-level annotation
  Output_source2 <- Output_source[, c("ASV", "Phylum")]
  Output_target2 <- Output_target[, c("ASV", "Phylum")]
  # Build summary table of Phylum–Phylum interactions
  summary <- data.frame(Source = Output_source2$Phylum,
                        Target = Output_target2$Phylum)
  # Generate adjacency matrix
  mig_data_filter <- as.matrix(as_adjacency_matrix(as_tbl_graph(summary)))
  
  # Remove overly broad or unwanted categories
  taxa_to_remove <- c("Kingdom:Bacteria")
  mig_data_filter2 <- mig_data_filter[!(rownames(mig_data_filter) %in% taxa_to_remove), ]
  mig_data_filter3 <- mig_data_filter2[, !(colnames(mig_data_filter2) %in% taxa_to_remove)]
  # Define Circos colors
  chord_color <- c(brewer.pal(12, "Set3"), brewer.pal(5, "Paired"))
  # Create Circos chord diagram
  chord <- chorddiag(
    data = mig_data_filter3,
    groupnamePadding = 25,
    groupPadding = 10,
    groupColors = chord_color,
    groupnameFontsize = 12,
    showTicks = FALSE,
    margin = 130,
    tooltipGroupConnector = "    &#x25B6;    ",
    chordedgeColor = "#B3B6B7",
    type = "directional"
  )
  # Save interactive and static visualizations
  filename_html <- paste0(treatment, ".html")
  saveWidget(chord, file = filename_html)
  webshot2::webshot(filename_html, file = filename_png, zoom = 6)
  webshot2::webshot(filename_html, file = filename_pdf, zoom = 6)
  # Store outputs for reuse
  Circos_list[[treatment]] <- list(
    chord = chord,
    mig_data_filter3 = mig_data_filter3,
    summary = summary
  )
}
