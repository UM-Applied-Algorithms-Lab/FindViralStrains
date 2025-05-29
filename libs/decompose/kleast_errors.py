import sys
import argparse
import gurobipy as gp
from gurobipy import GRB
from gurobipy import quicksum
import math
import networkx as nx
import itertools as it
import numpy as np
import flowpaths as fp
import time
import os
import matplotlib.pyplot as plt


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
                    flow = float(parts[-1]) 

                    # add edge 
                    if flow >= min_edge_weight:
                        graph.add_edge(u, v, flow=flow)
                except ValueError:
                    continue
                
    return graph  

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="""
        Decompose a network flow into a minimum number of weighted paths.
        This script uses the Gurobi ILP solver.
        """,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-t', '--threads', type=int, default=4, 
                       help='Number of threads to use for the Gurobi solver; use 0 for all threads (default 0).')
    parser.add_argument('-l', '--timelimit', type=int, default=100, 
                       help='time limit for Gurobi to solve one instance. (default 100)')
    parser.add_argument('-m', '--minpaths', type=int, default=1,
                       help='minimum number of paths to try (default 1).')
    parser.add_argument('-mc', '--mincount', type=int, default=0, 
                       help='minimum valid count on an edge (default 0)')
    parser.add_argument('-v', '--visualize', type=bool, default=False,
                       help='visualize the graph with matplotlib (default False).')
    
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
    requiredNamed.add_argument('-o', '--output', type=str, help='Output base filename', required=True)
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

def save_paths_to_file(paths, output_path, num_paths, runtime, mip_gap, objective_value):
    """Save path information to a text file in the specified format."""
    # Calculate total flow through all paths
    total_flow = sum(paths['weights'])
    
    with open(output_path, 'w') as f:
        f.write(f"Decomposition into {num_paths} paths\n")
        f.write(f"Runtime: {runtime:.2f} seconds\n")
        f.write(f"MIP Gap: {mip_gap:.6f}\n")
        f.write(f"Objective Value: {objective_value:.6f}\n")
        f.write(f"Number of Paths: {num_paths}\n")
        f.write("Paths and Weights:\n")
        
        for index, path in enumerate(paths['paths']):
            path_weight = paths['weights'][index]
            # Calculate fraction of total flow
            if total_flow > 0:
                fraction = path_weight / total_flow
            else:
                fraction = 0.0
            
            # Format the path with fraction first, then nodes
            path_str = " ".join(path)
            f.write(f"{fraction:.6f} {path_str}\n")
    
    print(f"INFO: Path details saved to {output_path}")

def save_graph_gml(graph, filepath):
    """Save graph in GML format with proper type handling"""
    try:
        
        nx.write_gml(graph, filepath, stringizer=str)
        print(f"Graph successfully saved to {filepath}")
    except Exception as e:
        print(f"Error saving GML: {str(e)}")
        raise



def save_graph_gexf(G, filepath):
    """Safely save graph as GEXF"""
    try:
        
        nx.write_gexf(
            G,
            filepath,
            encoding='utf-8',
            prettyprint=True,
            version='1.2draft'
        )
        print(f"Graph saved to {filepath}")
    except Exception as e:
        print(f"Error saving GEXF: {str(e)}")
        raise
def draw_labeled_multigraph(G, attr_name, ax=None, decimal_places= 2, paths = None):
    """
    Draw a multigraph with edge labels showing flow values.
    
    Parameters:
    - G: NetworkX graph
    - attr_name: Edge attribute to display
    - ax: Matplotlib axis (optional)
    - decimal_places: Number of decimal places to round to (default: 1)
    """
    # Connection styles for curved edges
    connectionstyle = [f"arc3,rad={r}" for r in it.accumulate([0.15] * 4)]

    # topologic sort the graph to ensure a proper layout
    list_of_nodes = list(nx.topological_sort(G))

    
    # Calculate dynamic figure size based on graph complexity
    num_nodes = graph.number_of_nodes()

    # Get graph layout
    pos = nx.nx_pydot.graphviz_layout(G, prog="sfdp")

    # place 0 and 1 at the furthest ends of the graph
    if '0' in pos and '1' in pos:
        pos['0'] = (min(pos.values(), key=lambda x: x[0])[0] - 1, pos['0'][1])
        pos['1'] = (max(pos.values(), key=lambda x: x[0])[0] + 1, pos['1'][1])
    
        # place both nodes at the same height, middway between the top and bottom of the graph
        min_y = min(y for _, y in pos.values())
        max_y = max(y for _, y in pos.values())
        mid_y = (min_y + max_y) / 2
        pos['0'] = (pos['0'][0], mid_y)
        pos['1'] = (pos['1'][0], mid_y)
    
    # place all other nodes in between the source and sink
    for node in list_of_nodes:
        if node != '0' and node != '1':
            # find the x position of the node by finding the average x position of its neighbors
            x_pos = np.mean([pos[neighbor][0] for neighbor in G.neighbors(node) if neighbor in pos])
            pos[node] = (x_pos, pos[node][1])  # Keep the y position the same
            
        
    
    # Adjust element sizes based on graph size
    font_size = max(8, 12 - math.log(num_nodes + 1))
    
    # Draw nodes and edges
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size= 300)
    nx.draw_networkx_labels(G, pos, font_size=font_size, ax=ax)
    nx.draw_networkx_edges(
                    G, pos, 
                    edge_color="grey", 
                    width=1.2,             # Thicker edges
                    connectionstyle=connectionstyle, 
                    ax=ax
                    )
    
    
    # Draw paths with thicker edges in red
    colors = ['r', 'b', 'g', 'c', 'm']                       # Colors for different paths
    for index, path in enumerate(paths['paths']):
        path_weight = paths['weights'][index]
        for i in range(len(path) - 1):
            u, v = path[i], path[i + 1]
            if graph.has_edge(u, v):
                
                
                nx.draw_networkx_edges(
                    G, pos, 
                    edgelist=[(u, v)], 
                    edge_color=colors[index % len(colors)],  # Cycle through colors
                    width=1.5,                               # Thicker edges for paths
                    connectionstyle=connectionstyle,
                    ax=ax
                )
            
                

    # Create edge labels with rounded values
    labels = {
        tuple(edge): f"{round(attrs[attr_name], decimal_places)}"
        for *edge, attrs in G.edges(keys=True, data=True)
    }
    
    # Draw edge labels
    nx.draw_networkx_edge_labels(
        G,
        pos,
        edge_labels=labels,
        font_size=9, 
        font_color="black",  
        font_weight="bold",         # Bold text
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "alpha": 0.7,            # Semi-transparent white background
            "edgecolor": "none"
        },
        rotate=False,   # Keep text horizontal
        label_pos=0.5,  # Center of edge
        connectionstyle=connectionstyle,
        ax=ax,
        horizontalalignment="center",  # Center text
        verticalalignment="center"     # Center text
    )



def visualize_and_save_graph(graph, output_path, num_paths, base_size=10, paths = None):
    """Visualize the graph and save to file."""
    # Calculate dynamic figure size based on graph complexity
    num_nodes = graph.number_of_nodes()
    num_edges = graph.number_of_edges()
    
    # Logarithmic scaling to handle large graphs
    scale_factor = math.log(num_nodes + num_edges + 1) * 0.75
    figsize = (base_size * scale_factor, base_size * scale_factor)
    

    fig, ax = plt.subplots(figsize = figsize)
    draw_labeled_multigraph(graph, 'flow', ax=ax, paths=paths)
    ax.set_title(f"Top {num_paths} Paths with Least Absolute Errors")
    
    visualization_file = f"{output_path}_visualization.pdf"
    plt.savefig(visualization_file, dpi=300, bbox_inches='tight')
    print(f"INFO: Visualization saved to {visualization_file}")



def generate_output_files(base_output_path, graph, max_paths, min_paths=1, visualize=False):
    """Generate output files for all path counts from max_paths down to min_paths."""
    # Extract the base filename without extension
    base_name = os.path.splitext(base_output_path)[0]

    
    for num_paths in range(max_paths, min_paths - 1, -1):
        # Create the specific output filename
        output_path = f"{base_name}_{num_paths}.paths"
        
        start_time = time.time()
        
        # Perform k-least errors analysis for current number of paths
        k_least = fp.kLeastAbsErrors(G=graph, k=num_paths, flow_attr='flow')
        k_least.solve()
        paths = k_least.get_solution(remove_empty_paths=True)
        
        # Get solver statistics
        runtime = time.time() - start_time
        mip_gap = k_least.model.MIPGap if hasattr(k_least, 'model') else 1.0
        objective_value = k_least.model.ObjVal if hasattr(k_least, 'model') else 0.0


        if visualize:
            # Visualize the graph
            visualize_and_save_graph(graph, output_path, num_paths, paths = paths)


        # Save path information
        save_paths_to_file(
            paths, 
            output_path, 
            num_paths,
            runtime,
            mip_gap,
            objective_value
        )

    
if __name__ == '__main__':
    # Parse command line arguments
    args = parse_arguments()
    print_thread_info(args.threads)

    # Read the input graph
    graph = read_graph_to_networkx(args.input, min_edge_weight=args.mincount)

    # Generate output files for all path counts from max_paths down to 1
    generate_output_files(args.output, graph, args.maxpaths, args.minpaths, visualize=args.visualize)

    print("INFO: Processing completed.")