import sys
import argparse
import gurobipy as gp
from gurobipy import GRB
from gurobipy import quicksum
import math
import networkx as nx
import itertools as it
import numpy as np
import matplotlib.pyplot as plt
import flowpaths as fp




def read_graph_to_networkx(file_path, min_edge_weight=0):
    """
    Reads a graph from a file and returns it as a NetworkX MultiDiGraph with flow attributes.
    """

    graph = nx.MultiDiGraph()
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            parts = line.split()
            if len(parts) >= 3:
                try:
                    u, v = parts[0], parts[1]
                    flow = float(parts[-2]) 

                    # add edge 
                    if flow >= min_edge_weight:
                        graph.add_edge(u, v, flow=flow)
                except ValueError:
                    continue
                
    return graph  



def draw_labeled_multigraph(G, attr_name, ax=None, decimal_places=1):
    """
    Draw a multigraph with edge labels showing flow values.
    
    Parameters:
    - G: NetworkX graph
    - attr_name: Edge attribute to display
    - ax: Matplotlib axis (optional)
    - decimal_places: Number of decimal places to round to (default: 1)
    """
    # Connection styles for curved edges
    connectionstyle = [f"arc3,rad={r}" for r in it.accumulate([0.40] * 5)]
    
    # Get graph layout
    pos = nx.nx_pydot.graphviz_layout(G)
    
    # Draw nodes and edges
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=500)
    nx.draw_networkx_labels(G, pos, font_size=10, ax=ax)
    nx.draw_networkx_edges(
        G, pos, 
        edge_color="grey", 
        width=1.5,  # Thicker edges
        connectionstyle=connectionstyle, 
        ax=ax
    )
    
    # Create edge labels with rounded values
    labels = {
        tuple(edge): f"{round(attrs[attr_name], decimal_places)}"
        for *edge, attrs in G.edges(keys=True, data=True)
    }
    
    # Draw edge labels with improved formatting
    nx.draw_networkx_edge_labels(
        G,
        pos,
        edge_labels=labels,
        font_size=9,  # Slightly larger font
        font_color="black",  # Higher contrast
        font_weight="bold",  # Bold text
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "alpha": 0.7,  # Semi-transparent white background
            "edgecolor": "none"
        },
        rotate=False,  # Keep text horizontal
        label_pos=0.5,  # Center of edge
        connectionstyle=connectionstyle,
        ax=ax,
        horizontalalignment="center",  # Center text
        verticalalignment="center"     # Center text
    )


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="""
        Decompose a network flow into a minimum number of weighted paths.
        This script uses the Gurobi ILP solver.
        """,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-t', '--threads', type=int, default=0, 
                       help='Number of threads to use for the Gurobi solver; use 0 for all threads (default 0).')
    parser.add_argument('-l', '--timelimit', type=int, default=100, 
                       help='time limit for Gurobi to solve one instance. (default 100)')
    parser.add_argument('-m', '--minpaths', type=int, default=1,
                       help='minimum number of paths to try (default 1).')
    parser.add_argument('-mc', '--mincount', type=int, default=0, 
                       help='minimum valid count on an edge (default 0)')
    
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
    requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)
    requiredNamed.add_argument('-M', '--maxpaths', type=int,
                             help='maximum number of paths to try.', required=True)
    
    return parser.parse_args()

def print_thread_info(threads):
    """Print information about thread usage."""
    if threads == 0:
        print("INFO: Using as many threads as possible (up to 32) for the Gurobi solver")
    else:
        print(f"INFO: Using {threads} threads for the Gurobi solver")

def create_k_least_graph(graph, paths):
    """Create a MultiDiGraph from the k-least errors paths."""
    k_least_graph = nx.MultiDiGraph()
    k_least_graph.add_nodes_from(graph.nodes())  # Ensure all nodes are included
    
    for index, path in enumerate(paths['paths']):
        path_weight = paths['weights'][index]
        for i in range(len(path) - 1):
            u, v = path[i], path[i + 1]
            if graph.has_edge(u, v):
                k_least_graph.add_edge(u, v, flow=path_weight)
    
    return k_least_graph

def visualize_and_save_graph(graph, output_path, num_paths):
    """Visualize the graph and save to file."""
    fig, ax = plt.subplots()
    draw_labeled_multigraph(graph, 'flow', ax=ax)
    ax.set_title(f"Top {num_paths} Paths with Least Absolute Errors")
    
    visualization_file = f"{output_path}_visualization.pdf"
    plt.savefig(visualization_file, dpi=300, bbox_inches='tight')
    print(f"INFO: Visualization saved to {visualization_file}")
    plt.show()

def save_paths_to_file(paths, output_path, num_paths):
    """Save path information to a text file."""
    with open(output_path, 'w') as f:
        f.write(f"Top {num_paths} Paths with Least Absolute Errors\n")
        f.write("="*50 + "\n\n")
        for index, path in enumerate(paths['paths']):
            path_weight = paths['weights'][index]
            f.write(f"Path {index + 1} (weight: {path_weight:.4f}):\n")
            f.write(" -> ".join(path) + "\n\n")
    print(f"INFO: Path details saved to {output_path}")


if __name__ == '__main__':
    
    # Parse command line arguments
    args = parse_arguments()
    print_thread_info(args.threads)

    # Read the input graph
    graph = read_graph_to_networkx(args.input, min_edge_weight=args.mincount)

    # Perform k-least errors analysis
    k_least = fp.kLeastAbsErrors(G=graph, k=args.maxpaths, flow_attr='flow')
    k_least.solve()
    paths = k_least.get_solution(remove_empty_paths=True)

    # Create and visualize graph
    k_least_graph = create_k_least_graph(graph, paths)
    visualize_and_save_graph(k_least_graph, args.output, args.maxpaths)

    # Save path information
    save_paths_to_file(paths, args.output, args.maxpaths)

    print("INFO: Processing completed.")


    

    






