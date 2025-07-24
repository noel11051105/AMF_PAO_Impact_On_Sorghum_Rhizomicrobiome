#NetLSD Distance
import igraph as ig
import numpy as np
import netlsd
# List of your GraphML files
graphml_files = [
    'AMF.graphml',
    'Control.graphml',
    'EBPR_AMF.graphml',
    'EBPR.graphml',
    'Vermont_AMF.graphml',
    'Vermont.graphml'
]
# Dictionary to store the NetLSD descriptors
descriptors = {}
for file in graphml_files:
    # Load the graph
    graph = ig.Graph.Read_GraphML(file)
    
    # Get the adjacency matrix
    adjacency_matrix = np.array(graph.get_adjacency().data)
    
    # Compute the NetLSD descriptor
    lsd_descriptor = netlsd.heat(adjacency_matrix)
    
    # Store the descriptor
    descriptors[file] = lsd_descriptor

# Print the descriptors
for file, descriptor in descriptors.items():
    print(f"{file}: {descriptor}")

# Compare the descriptors of all pairs of graphs and calculate the Euclidean distance
for i in range(len(graphml_files)):
    for j in range(i + 1, len(graphml_files)):
        descriptor1 = descriptors[graphml_files[i]]
        descriptor2 = descriptors[graphml_files[j]]
        
        # Calculate the Euclidean distance
        distance = np.linalg.norm(descriptor1 - descriptor2)
        print(f"Distance between {graphml_files[i]} and {graphml_files[j]} descriptors: {distance}")