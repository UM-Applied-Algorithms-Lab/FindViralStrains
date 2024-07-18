import sys
import argparse
import gurobipy as gp
from gurobipy import GRB
import time
import os

def get_extremity(neighbors, extremity_type):
    extremity = None
    for node, node_neighbors in neighbors.items():
        if len(node_neighbors) == 0:
            print(node)
            if extremity != None:
                sys.exit(f"ERROR: input graph has more than one {extremity_type}")
            else:
                extremity = node
    if extremity == None:
        sys.exit(f"ERROR: The input graph has no {extremity_type} (and thus it is not a DAG)")
    return extremity

def read_input_counts(graphfile, mincount):
    lines = open(graphfile ,'r').readlines()
    graphList = []
    for i in range(0, len(lines)):
        if lines[i][0] == '#':
            graphList.append(i)
    graphList.append(len(lines))
    
    finalGraphs = []
    for i in range(0, len(graphList) - 1):
        num_nodes = 0
        edges = dict()
        count = dict()
        out_neighbors = dict()
        in_neighbors = dict()
        max_count = 0
        
        for y in range(graphList[i], graphList[i+1]):
            line = lines[y].strip()
            if line[0] == '#' or line == '':
                continue
            elements = line.split('\t')
            if len(elements) == 1:
                num_nodes = int(elements[0])
            elif len(elements) == 3:
                count_value = int(elements[2])
                count[(elements[0], elements[1])] = 0 if count_value < mincount else count_value
                max_count = max(max_count, count_value)
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
        
        if num_nodes != len(in_neighbors):
            sys.exit(f"ERROR: expecting {num_nodes} nodes, the input graph has {len(in_neighbors)} nodes")
        
        source = get_extremity(in_neighbors, 'source')
        sink = get_extremity(out_neighbors, 'sink')
        finalGraphs.append({
            'vertices': list(in_neighbors.keys()),
            'count': count,
            'out_neighbors': out_neighbors,
            'in_neighbors': in_neighbors,
            'source': source,
            'sink': sink,
            'max_count': max_count
        })
    
    return finalGraphs

def decompose_flow(vertices, count, out_neighbors, in_neighbors, source, sink, max_count, K, threads, timelimit, outputfilename):
    V = vertices
    E = set(count.keys())
    data = dict()
    W = 1
    
    try:
        T = [(i, j, k) for (i, j) in E for k in range(0, K)]
        SC = [k for k in range(0, K)]
        CE = [(i, j) for (i, j) in E if count[i, j] > 0]
        
        model = gp.Model("MFD")
        model.Params.LogToConsole = 0
        model.Params.Threads = threads
        print(f"INFO: Trying to decompose into {K} paths...")
        
        x = model.addVars(T, vtype=GRB.BINARY, name="x")
        w = model.addVars(SC, vtype=GRB.CONTINUOUS, name="w", lb=0)
        z = model.addVars(T, vtype=GRB.CONTINUOUS, name="z", lb=0)
        r = model.addVars(CE, vtype=GRB.CONTINUOUS, name="r", lb=0)
        eps = model.addVars(CE, vtype=GRB.CONTINUOUS, name="eps", lb=0)
        
        for k in range(0, K):
            for i in V:
                if i == source:
                    model.addConstr(sum(x[i, j, k] for j in out_neighbors[i]) == 1)
                elif i == sink:
                    model.addConstr(sum(x[j, i, k] for j in in_neighbors[i]) == 1)
                else:
                    model.addConstr(
                        sum(x[i, j, k] for j in out_neighbors[i]) - sum(x[j, i, k] for j in in_neighbors[i]) == 0
                    )
        
        for (i, j) in E:
            for k in range(0, K):
                model.addConstr(z[i, j, k] <= W * x[i, j, k])
                model.addConstr(w[k] - (1 - x[i, j, k]) * W <= z[i, j, k])
                model.addConstr(z[i, j, k] <= w[k])
        
        model.addConstr(sum(w[k] for k in range(0, K)) == 1.0)
        
        for k in range(0, K - 1):
            model.addConstr(w[k] >= w[k + 1])
        
        for (i, j) in CE:
            model.addConstr(r[i, j] - (1.0/count[i, j]) * eps[i, j] <= (1.0/count[i, j]) * sum(z[i, j, k] for k in range(0, K)))
            model.addConstr(r[i, j] + (1.0/count[i, j]) * eps[i, j] >= (1.0/count[i, j]) * sum(z[i, j, k] for k in range(0, K)))
        
        for (i, j) in CE:
            for (k, l) in CE:
                if (i == k or j == l):
                    model.addConstr(r[i, j] == r[k, l])
        
        model.setObjective(sum(eps[i, j] for (i, j) in CE), GRB.MINIMIZE)
        model.Params.TimeLimit = timelimit
        model.optimize()
        
        print("Final MIP gap value: %f" % model.MIPGap)
        print('Obj: %g' % model.ObjVal)
        
        w_sol = [0] * len(range(0, K))
        x_sol = {}
        print("Model status: ", model.status)
        
        if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
            data['message'] = 'solved'
            data['runtime'] = model.Runtime
            data['mipgap'] = model.MIPGap
            data['objval'] = model.ObjVal
            
            for v in model.getVars():
                if 'w' in v.VarName:
                    elements = v.VarName.replace('[', ',').replace(']', ',').split(',')
                    k = int(elements[1])
                    w_sol[k] = v.x
                if 'x' in v.VarName:
                    elements = v.VarName.replace('[', ',').replace(']', ',').split(',')
                    i = elements[1]
                    j = elements[2]
                    k = int(elements[3])
                    x_sol[i, j, k] = v.x
            
            paths = extract_paths(x_sol, source, sink, out_neighbors, K)
            
            output_filename = f"{outputfilename}_{K}_paths.txt"
            with open(output_filename, 'w') as outputfile:
                outputfile.write(f"Decomposition into {K} paths\n")
                outputfile.write(f"Runtime: {data['runtime']:.2f} seconds\n")
                outputfile.write(f"MIP Gap: {data['mipgap']:.6f}\n")
                outputfile.write(f"Objective Value: {data['objval']:.6f}\n")
                outputfile.write(f"Number of Paths: {len(paths)}\n")
                outputfile.write("Paths and Weights:\n")
                
                if paths is None:
                    outputfile.write("ERROR: Did not find any decomposition\n")
                else:
                    for k in range(len(w_sol)):
                        outputfile.write(f"{w_sol[k]:.6f} {' '.join(map(str, paths[k]))}\n")
            
            print(f"INFO: Data for {K} paths decomposition has been written to {output_filename}")
            
            data['weights'] = w_sol
            data['paths'] = paths
        
        if model.status == GRB.INFEASIBLE:
            data['message'] = 'unsolved'
    
    except (gp.GurobiError, AttributeError, Exception) as e:
        print('Caught: ' + str(e))
        sys.exit('Exiting.')
    
    print(data['message'])
    return data

def extract_paths(x, source, sink, out_neighbors, K):
    paths = []
    for k in range(0, K):
        vertex = source
        path = [vertex]
        while vertex != sink:
            for out_neighbor in out_neighbors[vertex]:
                if x[vertex, out_neighbor, k] == 1:
                    vertex = out_neighbor
                    break
            path.append(vertex)
        paths.append(path)
    return paths

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
    parser.add_argument('-l', '--timelimit', type=int, default=100, help='time limit for Gurobi to solve one instance.')
    parser.add_argument('-mc', '--mincount', type=int, default=50, help='minimum valid count on an edge')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
    requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)

    args = parser.parse_args()

    threads = args.threads
    if threads == 0:
        threads = os.cpu_count()//2
    print(f"INFO: Using {threads} threads for the Gurobi solver")

    # Set minK and maxK to fixed values
    minK = 1
    maxK = 3
    timelimit = args.timelimit

    graphList = read_input_counts(args.input, args.mincount)
    curr = 0

    for graph in graphList:
        count = graph['count']
        vertices = graph['vertices']
        in_neighbors = graph['in_neighbors']
        out_neighbors = graph['out_neighbors']
        source = graph['source']
        sink = graph['sink']
        max_count = graph['max_count']

        for K in range(minK, maxK + 1):
            start = time.time()
            data = decompose_flow(vertices, count, out_neighbors, in_neighbors, source, sink, max_count, K, threads, timelimit, args.output)
            if data['message'] == "solved":
                end = time.time()
                print(f"INFO: Found a decomposition into {K} paths")
                # The output is now written within the decompose_flow function
                # Stop after finding a valid decomposition
        curr += 1

    print("INFO: Processing completed.")
