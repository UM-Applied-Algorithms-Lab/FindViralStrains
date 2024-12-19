import sys
import networkx as nx

def calculate_single_shortest_path(file_path, source, sink):
    graph = nx.DiGraph()

    # Parse the file
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Parse the number of nodes from the first non-comment line
    for line in lines:
        line = line.strip()
        if line and not line.startswith('#'):
            num_nodes = int(line)
            break

    # Parse edges
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue  # Skip empty lines and comments
        parts = line.split()
        if len(parts) == 3:
            u, v, weight = parts
            graph.add_edge(u, v, weight=int(weight))

    try:
        # Calculate the shortest path between the given source and sink
        length = nx.shortest_path_length(graph, source, sink, weight='weight')
        return length
    except nx.NetworkXNoPath:
        print(f"No path exists between {source} and {sink}")
        return None

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <graph_file> <source_node> <sink_node>")
        sys.exit(1)

    file_path = sys.argv[1]
    source_node = sys.argv[2]
    sink_node = sys.argv[3]
    
    try:
        shortest_path_length = calculate_single_shortest_path(file_path, source_node, sink_node)
        if shortest_path_length is not None:
            print(f"Shortest path length between {source_node} and {sink_node}: {shortest_path_length}")
    except Exception as e:
        print(f"An error occurred: {e}")

