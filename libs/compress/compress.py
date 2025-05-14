from collections import defaultdict, deque

class Graph:
    def __init__(self):
        self.nodes = set()
        self.edges = defaultdict(list)  # Adjacency list
        self.edge_data = {}  # Stores weight and sequence for each edge
    
    def add_edge(self, from_node, to_node, weight, sequence):
        self.nodes.add(from_node)
        self.nodes.add(to_node)
        self.edges[from_node].append(to_node)
        self.edge_data[(from_node, to_node)] = {
            'weight': weight,
            'sequence': sequence
        }
    
    def get_neighbors(self, node):
        return self.edges.get(node, [])
    
    def get_edge_data(self, from_node, to_node):
        return self.edge_data.get((from_node, to_node))
    
    def bfs(self, start_node):
        """Breadth-First Search traversal"""
        visited = set()
        queue = deque([start_node])
        traversal_order = []
        
        while queue:
            node = queue.popleft()
            if node not in visited:
                visited.add(node)
                traversal_order.append(node)
                for neighbor in self.get_neighbors(node):
                    if neighbor not in visited:
                        queue.append(neighbor)
        return traversal_order
    
    def dfs(self, start_node):
        """Depth-First Search traversal"""
        visited = set()
        stack = [start_node]
        traversal_order = []
        
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                traversal_order.append(node)
                # Push neighbors in reverse order to visit them in order
                for neighbor in reversed(self.get_neighbors(node)):
                    if neighbor not in visited:
                        stack.append(neighbor)
        return traversal_order

def load_graph_from_file(filename):
    graph = Graph()
    
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split(' ')               # Split by space for super source 
            parts_2 = line.strip().split('\t')            # Split by tab for all inner nodes
            if len(parts_2) >= 4:      
                from_node = int(parts_2[0])
                to_node = int(parts_2[1])
                weight = int(parts_2[2])
                sequence = parts_2[3]
                graph.add_edge(from_node, to_node, weight, sequence)    
            elif len(parts) == 3:                                       # length of super source is 3
                from_node = int(parts[0])
                to_node = int(parts[1])
                weight = int(parts[2])
                sequence = None          
                graph.add_edge(from_node, to_node, weight, sequence) . 
               
    
    return graph

if __name__ == "__main__":
    graph = load_graph_from_file('/Users/joserodriguez/54840891_S15_L001.super.wg')
    
    for int in range(0, 400):
       print(f"\nNeighbors of node {int}:", graph.get_neighbors(int))
   