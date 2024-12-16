from collections import defaultdict

def analyze_directed_graph(file_path):
    # Dictionaries to store in-degree and out-degree counts
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    
    # Read the file and process edges
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("#") or not line:  # Skip comments and empty lines
                continue

            # Split the line into components
            parts = line.split()
            if len(parts) < 3:
                continue  # Skip lines that don't have enough columns

            # Extract nodes and directions
            node1, node2 = parts[0], parts[1]
            node1_base = node1.rstrip('+-')  # Remove "+" or "-" for the base node ID
            node2_base = node2.rstrip('+-')
            node1_direction = node1[-1]  # "+" or "-"
            node2_direction = node2[-1]

            # Determine the directed relationship
            # If node1 ends with "+" and node2 ends with "-", treat it as node1 -> node2
            if node1_direction == "+" and node2_direction == "-":
                out_degree[node1_base] += 1
                in_degree[node2_base] += 1
            elif node1_direction == "-" and node2_direction == "+":
                out_degree[node2_base] += 1
                in_degree[node1_base] += 1
            else:
                # If directions are the same, treat as undirected edge for safety
                out_degree[node1_base] += 1
                out_degree[node2_base] += 1
                in_degree[node1_base] += 1
                in_degree[node2_base] += 1

    # Calculate statistics
    total_nodes = len(set(in_degree.keys()).union(out_degree.keys()))
    single_path_nodes = [
        node for node in set(in_degree.keys()).union(out_degree.keys())
        if in_degree[node] == 1 or out_degree[node] == 1
    ]

    return total_nodes, single_path_nodes

# Example usage
file_path = "/home/mikhail/Code/MFD-ILP/FindViralStrains/output/NoRefTest/mg/E1250_S84_L001/out.mg_subgraphs/graph_0.mg"
total_nodes, single_path_nodes = analyze_directed_graph(file_path)
print(f"Total number of nodes: {total_nodes}")
print(f"Number of nodes with only one path in or out: {len(single_path_nodes)}")
