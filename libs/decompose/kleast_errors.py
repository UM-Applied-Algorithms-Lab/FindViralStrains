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

def save_paths_to_file(paths, output_path, num_paths, runtime, mip_gap, objective_value, multigraph_decomposer=None):
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
            path_weight = paths['weights'][index] if total_flow > 0 else 0

            
            # Format the path appropriately based on whether we're dealing with a multigraph
            if multigraph_decomposer is not None:
                # For multigraphs, the path is a list of (u, v, key) tuples
                path_str = " ".join(f"{u}-{v}({key})" for u, v, key in path)
            else:
                # For regular graphs, the path is a list of nodes
                path_str = " ".join(path)
            
            f.write(f"{path_weight:.6f} {path_str}\n")
    
    print(f"INFO: Path details saved to {output_path}")


def draw_labeled_multigraph(G, attr_name, ax=None, decimal_places=2, paths=None):
    """
    Draw a multigraph with edge labels showing flow values and support for multiedge curvature.
    """
    if ax is None:
        ax = plt.gca()

    # Layout
    pos = nx.nx_pydot.graphviz_layout(G, prog='sfdp')

    # Manually adjust node '0' and '1' if present
    if '0' in pos and '1' in pos:
        # Original horizontal bounds
        x_coords = [x for x, _ in pos.values()]
        min_x, max_x = min(x_coords), max(x_coords)
        mid_y = sum([y for _, y in pos.values()]) / len(pos)

        # Set target spacing for '0' and '1'
        left_x = min_x - 1
        right_x = max_x + 1
        pos['0'] = (left_x, mid_y)
        pos['1'] = (right_x, mid_y)

        # Rescale other x-positions to lie strictly between '0' and '1'
        for node in pos:
            if node not in ['0', '1']:
                x, y = pos[node]
                # Scale x to [left_x + margin, right_x - margin]
                margin = 100
                scaled_x = ((x - min_x) / (max_x - min_x)) * (right_x - left_x - 2 * margin) + (left_x + margin)
                pos[node] = (scaled_x, y)

    # Font size based on node count
    num_nodes = G.number_of_nodes()
    font_size = max(8, 12 - math.log(num_nodes + 1))

    # Draw nodes and labels
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=300)
    nx.draw_networkx_labels(G, pos, font_size=font_size, ax=ax)

    # === Draw grey base edges with curvature ===
    for u, v, k in G.edges(keys=True):
        rad = 0.15 * (k if k % 2 == 1 else -k)  # Alternate direction
        style = f"arc3,rad={rad}"
        nx.draw_networkx_edges(
            G, pos,
            edgelist=[(u, v)],
            ax=ax,
            edge_color="grey",
            width=1.2,
            connectionstyle=style
        )

    # === Highlight paths if provided ===
    if paths and 'paths' in paths:
        colors = [
            "red", "blue", "green", "purple", "orange",
            "cyan", "magenta", "lime", "brown", "pink",
            "teal", "navy", "olive", "maroon", "darkgreen",
            "darkviolet", "gold"
        ]
        line_styles = ['solid', 'dashed', 'dotted', 'dashdot']
        line_widths = [5.00, 4.00, 3.00, 2.00]

        for idx, path in enumerate(paths['paths']):
            if all(isinstance(p, tuple) and len(p) == 3 for p in path):
                edge_list = path
            else:
                # Assume it's a list of nodes, default edge key 0
                edge_list = [(path[i], path[i+1], 0) for i in range(len(path)-1)]

            for u, v, k in edge_list:
                rad = 0.15 * (k if k % 2 == 1 else -k)
                style = f"arc3,rad={rad}"

                nx.draw_networkx_edges(
                    G, pos,
                    edgelist=[(u, v)],
                    edge_color=colors[idx % len(colors)],
                    width=line_widths[idx % len(line_widths)],
                    style=line_styles[idx % len(line_styles)],
                    connectionstyle=style,
                    ax=ax,
                    alpha=0.8,
                    arrows=True,
                    node_size=300,
                    arrowstyle='-|>',
                    arrowsize=15
                )

    # === Edge Labels ===
    edge_labels = {}
    for u, v, k, data in G.edges(keys=True, data=True):
        if attr_name in data:
            edge_labels[(u, v, k)] = f"{data[attr_name]:.{decimal_places}f}"

    # Place each label individually with correct curvature
    for (u, v, k), label in edge_labels.items():
        rad = 0.15 * (k if k % 2 == 1 else -k)
        try:
            nx.draw_networkx_edge_labels(
                G, pos,
                edge_labels={(u, v): label},
                font_size=9,
                font_color="black",
                font_weight="bold",
                bbox={
                    "boxstyle": "round",
                    "facecolor": "white",
                    "alpha": 0.7,
                    "edgecolor": "none"
                },
                rotate=False,
                label_pos=0.5,
                ax=ax,
                horizontalalignment="center",
                verticalalignment="center",
                connectionstyle=f"arc3,rad={rad}"
            )
        except Exception as e:
            print(f"Label error: {e}")

    ax.autoscale_view()
    plt.tight_layout()




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


def get_all_edges_for_node(graph, node):
    """
    Get all edges connected to a node in a NetworkX graph.
    Returns both incoming and outgoing edges.

    """
    # Convert node to string if it isn't already (since your keys are strings)
    node = str(node)
    outgoing = list(graph.adj.get(node, {}).keys())
    incoming = [u for u in graph.adj if node in graph.adj[u]]
    
    # Combine and format as edge tuples
    edges = [(node, v) for v in outgoing] + [(u, node) for u in incoming]
    
    return edges


def generate_output_files(base_output_path, graph, max_paths, min_paths=1, visualize=False):
    """Generate output files for all path counts from max_paths down to min_paths."""
    # Extract the base filename without extension
    base_name = os.path.splitext(base_output_path)[0]

    
    for num_paths in range(max_paths, min_paths - 1, -1):
        # Create the specific output filename
        output_path = f"{base_name}_{num_paths}.paths"
        
        start_time = time.time()

        # Get all edges to ignore for the current number of edges for source and sink nodes
        edges_to_ignore = get_all_edges_for_node(graph, "0") + get_all_edges_for_node(graph, "1")
    
        # Perform k-least errors analysis for current number of paths
        k_least = fp.kLeastAbsErrors(G=graph, k=num_paths, flow_attr='flow', elements_to_ignore=edges_to_ignore)
        k_least.solve()
        paths = k_least.get_solution(remove_empty_paths=True)
        

        # Get solver statistics
        runtime = time.time() - start_time
        mip_gap = k_least.model.MIPGap if hasattr(k_least, 'model') else 1.0
        objective_value = k_least.model.ObjVal if hasattr(k_least, 'model') else 0.0


        if visualize:
            # Visualize the graph
            visualize_and_save_graph(graph, output_path, num_paths, paths = paths)


        # see the type of graph
        decomposer = isinstance(graph, nx.MultiDiGraph)


        # Save path information
        save_paths_to_file(
            paths, 
            output_path, 
            num_paths,
            runtime,
            mip_gap,
            objective_value,
            multigraph_decomposer=decomposer
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