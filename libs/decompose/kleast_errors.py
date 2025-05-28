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

def generate_output_files(base_output_path, graph, max_paths, min_paths=1):
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

        # Create graph (without visualization)
        k_least_graph = create_k_least_graph(graph, paths)

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
    generate_output_files(args.output, graph, args.maxpaths, args.minpaths)

    print("INFO: Processing completed.")