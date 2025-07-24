# AMF_PAO_Impact_On_Sorghum_Rhizomicrobiome
Synergistic enhancement of Sorghum bicolor nutrient uptake and growth by EBPR microbiomes and AM fungi

This repository includes RData files and associated outputs used for univariate statistical analyses and 16S rRNA-based microbiome profiling, focused on plant biometrics, soil chemistry, and rhizosphere microbial communities.

Figure 1
- Fig1_codes.R: Contains statistical testing workflows and graphical outputs based on datasets related to plant biometrics (shoot height, dry mass, stalk diameter, root length, number of adventitious roots, AMF colonization rate), soil properties (moisture content, total and easily extractable glomalin-related soil proteins [GRSP], total carbon and nitrogen, organic matter), and elemental composition in sorghum roots (Al, Ca, Fe, K, Mg, P, S, Zn).
- Fig1.RData: Stores raw data and associated analysis outputs presented in Figure 1.
  
Figure 3
- Fig3_codes.R: Includes outputs from statistical testing and 16S microbiome analyses, such as relative abundance profiles of rhizomicrobial phyla and genera, phylum-level comparisons, NMDS ordination for treatment group clustering, and alpha diversity metrics (richness, evenness, Shannon index).
- Fig3.RData: Contains data objects and visual results associated with Figure 3.
  
Figure 4
- Fig4_codes.R: Provides results from differential abundance analysis using ALDEx2, Spearman correlation between microbial taxa and plant nutrient uptake, and microbial associations with alpha diversity indices. Also includes PICRUSt2 outputs for functional trait inference in rhizosphere microbial communities.
- Fig4.RData: Stores raw datasets and processed outputs for Figure 4.
- KO.tsv: Tab-delimited file containing KEGG Orthology (KO) abundance data.
  
Figure 5
- Fig5_codes.R: Contains results from microbial network analyses, including files for Cytoscape visualization and network structure comparisons using Portrait Divergence and NetLSD metrics.
- PortraitDivergence.py: Python script for executing NetLSD-based network similarity assessments.
- NetNSL.py: Python codes required to run NetNSL analyses.
- Fig5.RData: Includes analysis results and visualization-ready data featured in Figure 5.
