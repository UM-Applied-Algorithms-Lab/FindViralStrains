import sys
import argparse
import gurobipy as gp
from gurobipy import GRB
import logging
import math

def read_input_counts(graph_file_src, min_edge_weight):
    """ Reads a weighted graph file and extracts relevant data. """
    with open(graph_file_src, "r") as graph_file:
        lines = graph_file.readlines()
    
    edge_weights = {}
    out_neighbors = {}
    in_neighbors = {}
    max_edge_weight = 0

    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        elements = line.split()
        if len(elements) == 3:
            u, v, weight = elements[0], elements[1], int(elements[2])
            edge_weights[(u, v)] = max(weight, min_edge_weight)
            max_edge_weight = max(max_edge_weight, weight)
            
            if u not in out_neighbors:
                out_neighbors[u] = []
            if v not in in_neighbors:
                in_neighbors[v] = []
            
            out_neighbors[u].append(v)
            in_neighbors[v].append(u)
    
    source = next(node for node in out_neighbors if node not in in_neighbors)
    sink = next(node for node in in_neighbors if node not in out_neighbors)
    
    return {
        'vertices': list(in_neighbors.keys()),
        'count': edge_weights,
        'out_neighbors': out_neighbors,
        'in_neighbors': in_neighbors,
        'source': source,
        'sink': sink,
        'max_count': max_edge_weight
    }

class Encode_Modified:
    def __init__(self, graph_data, num_paths, epsilon, timeout, threads):
        self.n = len(graph_data['vertices'])
        self.E = list(graph_data['count'].keys())
        self.F = graph_data['count']
        self.source = graph_data['source']
        self.target = graph_data['sink']
        self.k = num_paths
        self.w_max = graph_data['max_count']
        self.epsilon = epsilon
        self.timeout = timeout
        self.threads = threads
        
        self.edge_vars = {}
        self.weights = {}
        self.error_vars = {}
        self.model = self.create_solver()
    
    def create_solver(self):
        env = gp.Env(empty=True)
        env.setParam('OutputFlag', 0)
        env.setParam('TimeLimit', self.timeout)
        env.setParam('Threads', self.threads)
        env.start()
        return gp.Model("Modified_ILP", env=env)
    
    def encode(self):
        edge_indexes = [(u, v, i) for i in range(self.k) for (u, v) in self.E]
        path_indexes = [i for i in range(self.k)]
        
        self.edge_vars = self.model.addVars(edge_indexes, vtype=GRB.BINARY, name='x')
        self.weights = self.model.addVars(path_indexes, vtype=GRB.CONTINUOUS, name='w', lb=0)
        self.error_vars = self.model.addVars(self.E, vtype=GRB.CONTINUOUS, name='eps', lb=0)
        
        for i in range(self.k):
            self.model.addConstr(self.edge_vars.sum(self.source, '*', i) == 1)
            self.model.addConstr(self.edge_vars.sum('*', self.target, i) == 1)
            
        for u, v in self.E:
            flow_sum = self.edge_vars.sum(u, v, '*')
            self.model.addConstr(self.F[(u, v)] - flow_sum <= self.error_vars[u, v])
            self.model.addConstr(flow_sum - self.F[(u, v)] <= self.error_vars[u, v])
        
        self.model.setObjective(self.error_vars.sum(), GRB.MINIMIZE)
    
    def solve(self):
        self.model.optimize()
        return self.model.status, self.extract_solution()
    
    def extract_solution(self):
        paths = []
        for i in range(self.k):
            path = []
            u = self.source
            while u != self.target:
                next_edges = [v for v in self.model.getVars() if v.X > 0.9 and v.VarName.startswith(f"x[{u},")]
                if not next_edges:
                    break
                v = next_edges[0].VarName.split(',')[1].strip(']')
                path.append(v)
                u = v
            paths.append(path)
        return paths

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Modified ILP Flow Decomposition")
    parser.add_argument('-i', '--input', required=True, help='Input graph file')
    parser.add_argument('-p', '--paths', type=int, required=True, help='Number of paths')
    parser.add_argument('-e', '--epsilon', type=float, default=0.25, help='Error tolerance')
    parser.add_argument('-t', '--timeout', type=int, default=300, help='Solver time limit')
    parser.add_argument('-th', '--threads', type=int, default=4, help='Number of solver threads')
    args = parser.parse_args()

    graph_data = read_input_counts(args.input, min_edge_weight=1)
    encoder = Encode_Modified(graph_data, args.paths, args.epsilon, args.timeout, args.threads)
    encoder.encode()
    status, solution = encoder.solve()

    if status == GRB.OPTIMAL:
        print("Optimal solution found:")
        for path in solution:
            print(" -> ".join(path))
    else:
        print("No optimal solution found within time limit.")
