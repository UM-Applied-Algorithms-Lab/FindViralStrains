import sys
import argparse
import gurobipy as gp
from gurobipy import GRB
from gurobipy import quicksum
import time
import os
import math


def get_extremity(neighbors, extremity_type):
    """TODO: deprecate this function, and use super-source and super-sink files 
    Returns the 'extremity', or the node with no neighbors in the given neighbor set
    Function will exit the program if multiple extremities are found. 
    For example, this function can identify sources or sinks, depending on the neighbors collection given.
    """
    
    #filter to find any nodes without neighbors 
    extremity_list = [node for (node, neighbor_nodes) in neighbors.items() if len(neighbor_nodes) == 0]
    
    # any number of extremities other than 1 is an error
    if len(extremity_list) == 0:
        sys.exit(f"ERROR: The input graph has no {extremity_type} (and thus it is not a DAG)")
    elif len(extremity_list) > 1:
        sys.exit(f"ERROR: input graph has more than one {extremity_type}")
    else:
        return extremity_list[0]
        
    # extremity = None
    # for node, node_neighbors in neighbors.items():
    #     if len(node_neighbors) == 0:
    #         print(node)
    #         if extremity != None:
    #             sys.exit(f"ERROR: input graph has more than one {extremity_type}")
    #         else:
    #             extremity = node
    # if extremity == None:
    #     sys.exit(f"ERROR: The input graph has no {extremity_type} (and thus it is not a DAG)")
    # return extremity
    
    
def read_source_sink_file(source_sink_file_src):
    """reads the super source or super sink file, and gets the name of the source/sink node
    """
    with open(source_sink_file_src, "r") as source_sink_file:
        return source_sink_file.readline()


def read_input_counts(graph_file_src, min_edge_weight):
    """ Read through the input weighted graph file, and returns the data as a dict containing:
    list of vertices,
    weights of the edges,
    out neighbors,
    in neighbors,
    the source node,
    the sink node,
    and the maximum weight along all edges
    
    """
    finalGraphs = []
    with open(graph_file_src, "r") as graph_file:
    
        lines = graph_file.readlines()
        graphList = []
        for i in range(0, len(lines)):
            if lines[i][0] == '#':
                graphList.append(i)
        graphList.append(len(lines))

        for i in range(0, len(graphList) - 1):
            num_edges = 0
            edge_weights = dict()
            out_neighbors = dict()
            in_neighbors = dict()
            max_edge_weight = 0

            for y in range(graphList[i], graphList[i+1]):
                line = lines[y].strip()
                if len(line) == 0 or line[0] == '#':
                    continue
                
                elements = line.split()
                if len(elements) == 1:
                    num_edges = int(elements[0])
                elif len(elements) == 3 or len(elements) == 6:  # accepts lines with 3 or 6 elements
                    edge_weight_value = float(elements[3])
                    edge_weights[(elements[0], elements[1])] = 0 if edge_weight_value < min_edge_weight else edge_weight_value    #this line turns all counts below min_count to zero! defaults to 0, TODO
                    max_edge_weight = max(max_edge_weight, edge_weight_value) #just used for a global max edges
                    if elements[0] not in out_neighbors:
                        out_neighbors[elements[0]] = []
                    out_neighbors[elements[0]].append(elements[1])
                    if elements[1] not in in_neighbors:
                        in_neighbors[elements[1]] = []
                    in_neighbors[elements[1]].append(elements[0])
                    if elements[0] not in in_neighbors:
                        in_neighbors[elements[0]] = []
                    if elements[1] not in out_neighbors:
                        out_neighbors[elements[1]] = []
                else:
                    sys.exit("ERROR: input file contains an ill-formatted line")

            if num_edges != len(in_neighbors) - 2:
                sys.exit(f"ERROR: expecting {num_edges} edges, the input graph has {len(in_neighbors)} edges")

            source = get_extremity(in_neighbors, 'source')
            sink = get_extremity(out_neighbors, 'sink')
            finalGraphs.append({
                'vertices': list(in_neighbors.keys()),
                'count': edge_weights,
                'out_neighbors': out_neighbors,
                'in_neighbors': in_neighbors,
                'source': source,
                'sink': sink,
                'max_count': max_edge_weight
            })

    return finalGraphs


def decompose_flow(vertices, count, out_neighbors, in_neighbors, source_node_name, sink_node_name, \
    max_count, num_paths, num_threads, time_limit, output_file_name):
    
    edges = list(set(count)) # Pull edge weights, TODO check if this is ok #
    out_neighbors = {k: list(set(v)) for k, v in out_neighbors.items()}
    in_neighbors = {k: list(set(v)) for k, v in in_neighbors.items()}
    output_data = dict()
    W = 1 # Max flow

    try:
        T = [(from_node, to_node, k) for (from_node, to_node) in edges for k in range(0, num_paths)]    #collection of (i)node->(j)node edges for each (k)path
        SC = list(range(0, num_paths))

        print(f"INFO: Setting up model...")
        model = gp.Model("MFD")
        model.Params.LogToConsole = 0
        model.Params.Threads = num_threads

#-------------------------------------------------------------------------------------------------
# VARS
#-------------------------------------------------------------------------------------------------
# lb is lower bound for a variable
        x = model.addVars(T, vtype=GRB.BINARY, name="x") # binary variable to show if path k uses edge ij
        w = model.addVars(SC, vtype=GRB.CONTINUOUS, name="w", lb=0) # Fractional flow of path k
        z = model.addVars(T, vtype=GRB.CONTINUOUS, name="z", lb=0) # Fractional flow of path k carried on edge
        f = model.addVars(edges, vtype=GRB.CONTINUOUS, name="f", lb=0) # Fit flow on edge ij
        
        # New variable for path flow error
        path_error = model.addVars(edges, num_paths, vtype=GRB.CONTINUOUS, name="path_error", lb=0)
        epsilon = model.addVars(edges, vtype=GRB.CONTINUOUS, name="eps", lb=0) # Min error
#-------------------------------------------------------------------------------------------------
# Constraints
#-------------------------------------------------------------------------------------------------
       
        print(f"Creating {len(vertices) * num_paths} flow conservation constraints")
        #constraints 7,8
        for path_idx in range(0, num_paths):
            for vertex in vertices:
                if vertex == source_node_name:
                    # the edge weights outgoing from a source must equal 1 (can't take multiple edges out of the source on a single path)
                    model.addConstr(sum(x[vertex, neighbor_node, path_idx] for neighbor_node in out_neighbors[vertex]) == 1)
                elif vertex == sink_node_name:
                    # the edge weights incoming to a sink node must equal 1 (can't take multiple edges into the sink on a single path)
                    model.addConstr(sum(x[neighbor_node, vertex, path_idx] for neighbor_node in in_neighbors[vertex]) == 1)
                else:
                    # the path coming into a node must come out of that node.
                    model.addConstr(
                        sum(x[vertex, neighbor_node, path_idx] for neighbor_node in out_neighbors[vertex]) == \
                            sum(x[neighbor_node, vertex, path_idx] for neighbor_node in in_neighbors[vertex])
                    )
        
        print(f"Creating {len(edges) * num_paths} linearization constraints")
        #constraints 9,10,11
        for (vertex_from, vertex_to) in edges:
            for path_idx in range(0, num_paths):
                model.addConstr(z[vertex_from, vertex_to, path_idx] <= W * x[vertex_from, vertex_to, path_idx])
                model.addConstr(w[path_idx] - (1 - x[vertex_from, vertex_to, path_idx]) * W <= z[vertex_from, vertex_to, path_idx])
                model.addConstr(z[vertex_from, vertex_to, path_idx] <= w[path_idx])


        for path_idx in range(0, num_paths - 1):
            model.addConstr(w[path_idx] >= w[path_idx + 1])

        print(f"Creating {2 * len(edges) * num_paths} path flow error constraints")
        
        
        #constraints 3,4,5,6,12,13
        for (vertex_from, vertex_to) in edges:
            
            total_outgoing_count = sum(count[vertex_from, neighbor] for neighbor in out_neighbors[vertex_from])
            total_incoming_count = sum(count[neighbor, vertex_to] for neighbor in in_neighbors[vertex_to])
            total_flow_in = sum(f[neighbor, vertex_to] for neighbor in in_neighbors[vertex_to])
            total_flow_out = sum(f[vertex_from, neighbor] for neighbor in out_neighbors[vertex_from])
            
            #constraint 5
            model.addConstr(total_incoming_count*f[vertex_from, vertex_to] - count[vertex_from, vertex_to]*total_flow_in  <= total_incoming_count * epsilon[vertex_from, vertex_to])
           
            #constraint 6
            model.addConstr(count[vertex_from, vertex_to]*total_flow_in - total_incoming_count*f[vertex_from, vertex_to]  <= total_incoming_count * epsilon[vertex_from, vertex_to])
            
            #constraint 3
            model.addConstr(total_outgoing_count*f[vertex_from, vertex_to] - count[vertex_from, vertex_to]*total_flow_out <= total_outgoing_count * epsilon[vertex_from, vertex_to])

            #constraint 4
            model.addConstr(count[vertex_from, vertex_to]*total_flow_out - total_outgoing_count*f[vertex_from, vertex_to]  <= total_outgoing_count * epsilon[vertex_from, vertex_to])
        

            #13 Actual flow - expected flow <= path flow error
            model.addConstr(
                    sum(z[vertex_from, vertex_to, path_idx] for path_idx in range (0, num_paths)) - f[vertex_from, vertex_to] <= path_error[vertex_from, vertex_to, path_idx],
                    name=f"path_error_pos_{vertex_from}_{vertex_to}_{path_idx}"
                )
            #12 Actual flow - expected flow <= path flow error
            model.addConstr(
                    f[vertex_from, vertex_to] - sum(z[vertex_from, vertex_to, path_idx] for path_idx in range (0, num_paths)) <= path_error[vertex_from, vertex_to, path_idx],
                    name=f"path_error_pos_{vertex_from}_{vertex_to}_{path_idx}"
                )

        print("Starting to create learned flow constraints")
        
        #constraint 2
        for vertex in vertices:
            if (vertex != "0") and (vertex != "1"):
                model.addConstr(sum(f[vertex, neighbor_1] for neighbor_1 in out_neighbors[vertex]) == 
                            sum(f[neighbor_2, vertex] for neighbor_2 in in_neighbors[vertex]) )
           
        #Constraint 1
        model.addConstr(sum(f["0", neighbor_1] for neighbor_1 in out_neighbors["0"]) == 1)
        model.addConstr(sum(f[neighbor_2, "1"] for neighbor_2 in in_neighbors["1"]) == 1)

        
        print(f"Added flow per read constraints")

#-------------------------------------------------------------------------------------------------
# OBJECTIVE
#-------------------------------------------------------------------------------------------------
        model.setObjective(
            quicksum(path_error[vertex_from, vertex_to, path_idx] 
                     for (vertex_from, vertex_to) in edges
                     for path_idx in range(num_paths))+ quicksum(epsilon[vertex_from, vertex_to] 
             for (vertex_from, vertex_to) in edges) , 
            GRB.MINIMIZE
        )

        
        model.Params.TimeLimit = time_limit
        print(model.Params.TimeLimit)
        print(f"INFO: Trying to decompose into {num_paths} paths...")
        model.optimize()

    # Debugging statements after optimization
        for (vertex_from, vertex_to) in edges:
            total_outgoing_count = sum(count[vertex_from, neighbor] for neighbor in out_neighbors[vertex_from])
            print([(vertex_from, neighbor )for neighbor in out_neighbors[vertex_from]])
            total_incoming_count = sum(count[neighbor, vertex_to] for neighbor in in_neighbors[vertex_to])
            total_flow_in = sum(f[neighbor, vertex_to].x for neighbor in in_neighbors[vertex_to])
            total_flow_out = sum(f[vertex_from, neighbor].x for neighbor in out_neighbors[vertex_from])

            # Debugging statements to check the values after optimization
            print(f"vertex_from: {vertex_from}, vertex_to: {vertex_to}")
            print(f"total_outgoing_count: {total_outgoing_count}, total_incoming_count: {total_incoming_count}")
            print(f"total_flow_in: {total_flow_in}, total_flow_out: {total_flow_out}")
            print(f"count[vertex_from, vertex_to]: {count[vertex_from, vertex_to]}")
            print(f"epsilon[vertex_from, vertex_to]: {epsilon[vertex_from, vertex_to].x}")

            # Print the errors for each edge
            print(f"Edge {vertex_from} to {vertex_to} has error {epsilon[vertex_from, vertex_to].x}")

            # Print the actual flow for each edge
            print(f"Edge {vertex_from} to {vertex_to} has actual flow {f[vertex_from, vertex_to].x}")

            # Print the expected flow for each edge
            expected_flow = sum(z[vertex_from, vertex_to, path_idx].x for path_idx in range(num_paths))
            print(f"Edge {vertex_from} to {vertex_to} has expected flow {expected_flow}")

            
        
        print("Final MIP gap value: %f" % model.MIPGap)
        print('Obj: %g' % model.ObjVal)

        w_sol = [0] * len(range(0, num_paths))
        x_sol = {}
        r_sol = {}
        eps_sol = {}
        path_flow = {}
        print("Model status: ", model.status)

        if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
            output_data['message'] = 'solved'
            output_data['runtime'] = model.Runtime
            output_data['mipgap'] = model.MIPGap
            output_data['objval'] = model.ObjVal

            for model_var in model.getVars():
                if "w" in model_var.VarName:
                    elements = model_var.VarName.replace('[', ',').replace(']', ',').split(',')
                    path_idx = int(elements[1])
                    w_sol[path_idx] = model_var.x
                    
                if 'x' in model_var.VarName:
                    elements = model_var.VarName.replace('[', ',').replace(']', ',').split(',')
                    vertex_from = elements[1]
                    vertex_to = elements[2]
                    path_idx = int(elements[3])
                    x_sol[vertex_from, vertex_to, path_idx] = model_var.x
            
                if 'f' in model_var.VarName:
                    elements = model_var.VarName.replace('[', ',').replace(']', ',').split(',')
                    vertex_from = elements[1]
                    vertex_to = elements[2]
                    r_sol[vertex_from, vertex_to] = model_var.x
                if 'eps' in model_var.VarName:
                    elements = model_var.VarName.replace('[', ',').replace(']', ',').split(',')
                    vertex_from = elements[1]
                    vertex_to = elements[2]
                    eps_sol[vertex_from, vertex_to] = model_var.x
                if 'path_error' in model_var.VarName:
                    elements = model_var.VarName.replace('[', ',').replace(']', ',').split(',')
                    vertex_from = elements[1]
                    vertex_to = elements[2]
                    path_flow[vertex_from, vertex_to] = model_var.x
                
          
        
            paths = extract_paths(x_sol, source_node_name, sink_node_name, out_neighbors, num_paths)
            output_data['weights'] = w_sol
            output_data['paths'] = paths
            output_data['r_sol'] = r_sol
            output_data['eps_sol'] = eps_sol
            output_data['path_flow'] = path_flow
        if model.status == GRB.INFEASIBLE:
            output_data['message'] = 'unsolved'

    except (gp.GurobiError, AttributeError, Exception) as e:
        print('Caught: ' + str(e))
        sys.exit('Exiting.')

    print(output_data['message'])
    return output_data

#def extract_paths(x, source, sink, out_neighbors, K):
#    paths = []
#    for k in range(0, K):
#        vertex = source
#        path = [vertex]
#        print(path)
#        while vertex != sink:
#            for out_neighbor in out_neighbors[vertex]:
#                if x[vertex, out_neighbor, k] == 1:
#                    vertex = out_neighbor
#                    break
#            path.append(vertex)
#        paths.append(path)
#    return paths

def extract_paths(x, source, sink, out_neighbors, K):
    paths = []
    for k in range(0, K):
        vertex = source
        path = [vertex]
        print(f"Starting path extraction for path {k} from source {source} to sink {sink}")

        while vertex != sink:
            found = False
            for out_neighbor in out_neighbors[vertex]:
                if (vertex, out_neighbor, k) in x and math.isclose(x[vertex,out_neighbor, k],1,rel_tol=1e-5 ):
                # Use math.isclose, because gurobi sometimes returns numbers like 1.00000007 #
                    vertex = out_neighbor
                    path.append(vertex)
                    found = True
                    break
            
            if not found:
                print(f"Warning: No valid outgoing edge found for vertex {vertex} in path {k}")
                # Later on change thus to throw a real error to warn the user #
                break

        paths.append(path)

    return paths


def write_results(data, outputfilename, num_paths):

    output_filename = f"{outputfilename[:-4]}_{num_paths}.paths"
    with open(output_filename, 'w') as outputfile:
        outputfile.write(f"Decomposition into {num_paths} paths\n")
        outputfile.write(f"Runtime: {data['runtime']:.2f} seconds\n")
        outputfile.write(f"MIP Gap: {data['mipgap']:.6f}\n")
        outputfile.write(f"Objective Value: {data['objval']:.6f}\n")
        outputfile.write(f"Number of Paths: {len(data['paths'])}\n")
        outputfile.write("Paths and Weights:\n")
        print(f"INFO: Data for {num_paths} paths decomposition has been written to {output_filename}")

        if data["paths"] is None:
            outputfile.write("ERROR: Did not find any decomposition\n")
        else:
            for k in range(len(data["weights"])):
                outputfile.write(f"{data['weights'][k]:.6f} {' '.join(map(str, data['paths'][k]))}\n")


def write_graph_results(data, outputfilename, K, count):
    output_filename = f"{outputfilename}_{K}_paths.txt"
    with open(output_filename, 'w') as outputfile:
        for edge in count:
            outputfile.write(f"{edge[0]}\t{edge[1]}\t")
            try:
                outputfile.write(f'{data["r_sol"][edge]:.4f},{data["path_flow"][edge]:.4f},{data["eps_sol"][edge]:.4f},{count[edge[0], edge[1]]}\n')
            except KeyError:
                outputfile.write(f'x,x,{count[edge[0], edge[1]]}\n')


def write_table_line(output_file_src, inputfilename, data):
    with open(output_file_src, 'a') as output_file:
        output_file.write(str(data["objval"]) + '\t')
        if data["paths"] == None:
            output_file.write("ERROR: Did not find any decomposition")
        else:
            for k in range(0, len(data["weights"])):
                output_file.write(f"{data['weights'][k]}\t")
    print("INFO: Data has been written to the output file.")  # Add this line


def write_reset(output_file_src, input_file_src):
    with open(output_file_src, "a") as output_file:
        output_file.write(input_file_src + "\t")


def write_return(output_file_src):
    with open(output_file_src, 'a') as output_file:
        output_file.write("\n")


if __name__ == '__main__':
    # Argument parser
    parser = argparse.ArgumentParser(
        description="""
        Decompose a network flow into a minimum number of weighted paths.
        This script uses the Gurobi ILP solver.
        """,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-t', '--threads', type=int, default=0, help='Number of threads to use for the Gurobi solver; use 0 for all threads (default 0).')
    parser.add_argument('-l', '--timelimit', type=int, default=100, help='time limit for Gurobi to solve one instance. (default 100)')
    parser.add_argument('-m', '--minpaths', type=int, default=1,
                        help='minimum number of paths to try (default 1).')
    parser.add_argument('-mc', '--mincount', type=int, default=0, help='minimum valid count on an edge (default 0)')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
    requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)
    requiredNamed.add_argument('-M', '--maxpaths', type=int,
                               help='maximum number of paths to try.', required=True)
    args = parser.parse_args()

    threads = args.threads
    if threads == 0:
        print("INFO: Using as many threads as possible (up to 32) for the Gurobi solver")
    else:
        print(f"INFO: Using {threads} threads for the Gurobi solver")

    # Set minK and maxK to fixed values
    min_num_paths = args.minpaths
    max_num_paths = args.maxpaths
    time_limit = args.timelimit

    graphList = read_input_counts(args.input, args.mincount)

    for graph in graphList:
        count = graph['count']
        vertices = graph['vertices']
        in_neighbors = graph['in_neighbors']
        out_neighbors = graph['out_neighbors']
        source = graph['source']
        sink = graph['sink']
        max_count = graph['max_count']
        write_reset(args.output, args.input)
        for K in range(min_num_paths, max_num_paths + 1):
            start = time.time()
            data = decompose_flow(vertices, count, out_neighbors, in_neighbors, source, sink, max_count, K, threads, time_limit, args.output)
            if data['message'] == "solved":
                write_table_line(args.output, args.input, data)
                end = time.time()
                write_results(data, args.output, K)
                write_graph_results(data, args.output + ".graph", K, count)
                print(f"INFO: Found a decomposition into {K} paths")
                # The output is now written within the decompose_flow function
                # Stop after finding a valid decomposition
        write_return(args.output)
    print("INFO: Processing completed.")

