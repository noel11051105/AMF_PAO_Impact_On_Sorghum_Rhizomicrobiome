#Portrait Divergence
import networkx as nx
import sys
sys.path.append('/network-portrait-divergence/')
from portrait_divergence import portrait_divergence

# List of GraphML files
graphml_files = [
    "AMF.graphml",
    "Control.graphml",
    "EBPR_AMF.graphml",
    "EBPR.graphml",
    "Vermont_AMF.graphml",
    "Vermont.graphml"
]

# Loop through each pair of files and calculate portrait divergence
for i in range(len(graphml_files)):
    for j in range(i + 1, len(graphml_files)):
        # Load the graphs
        Ger1 = nx.read_graphml(graphml_files[i])
        Ger2 = nx.read_graphml(graphml_files[j])
        
        # Calculate and print the portrait divergence
        print(f"Djs({graphml_files[i]}, {graphml_files[j]}) =", portrait_divergence(Ger1, Ger2))