from collections import defaultdict, namedtuple
import sys

Edge = namedtuple('Edge', ['to', 'weight', 'seq'])

def read_graph(filename):
    forward_edges = defaultdict(list)  # node -> [Edge]
    reverse_edges = defaultdict(list)  # node -> [from_nodes]
    node_seqs = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) < 4:
                continue
                
            try:
                from_node = int(parts[0])
                to_node = int(parts[1])
                weight = int(parts[2])
                seq = ' '.join(parts[3:])
                
                edge = Edge(to_node, weight, seq)
                forward_edges[from_node].append(edge)
                reverse_edges[to_node].append(from_node)
                node_seqs[from_node] = seq
                if to_node not in node_seqs:
                    node_seqs[to_node] = ""
            except ValueError:
                continue
    
    return forward_edges, reverse_edges, node_seqs

def find_merge_candidates(forward_edges, reverse_edges):
    merge_candidates = []
    all_nodes = set(forward_edges.keys()).union(
        edge.to for edges in forward_edges.values() for edge in edges)
    
    for node in sorted(all_nodes):
        incoming_count = len(reverse_edges.get(node, []))
        outgoing_count = len(forward_edges.get(node, []))
        
        if incoming_count == 1 and outgoing_count == 1:
            source = reverse_edges[node][0]
            target = forward_edges[node][0].to
            merge_candidates.append((source, node, target))
    
    return merge_candidates

def merge_nodes(forward_edges, reverse_edges, node_seqs):
    changed = True
    while changed:
        changed = False
        candidates = find_merge_candidates(forward_edges, reverse_edges)
        
        for source, node, target in candidates:
            # Verify the nodes still meet the merge conditions (might have changed)
            if (len(reverse_edges.get(node, [])) == 1 and 
                len(forward_edges.get(node, [])) == 1 and
                reverse_edges[node][0] == source and
                forward_edges[node][0].to == target):
                
                # Get the edge weights
                edge1_weight = next(e.weight for e in forward_edges[source] if e.to == node)
                edge2_weight = forward_edges[node][0].weight
                
                # Calculate new average weight
                new_weight = (edge1_weight + edge2_weight) // 2
                
                # Combine sequences
                new_seq = node_seqs[source] + node_seqs[node][-1:]  # sequences overlap by k-1
                print(node_seqs[node][-1:])
                
                # Remove old edges
                forward_edges[source] = [e for e in forward_edges[source] if e.to != node]
                forward_edges[node] = []
                reverse_edges[node] = []
                reverse_edges[target] = [n for n in reverse_edges[target] if n != node]
                
                # Add new edge
                new_edge = Edge(target, new_weight, new_seq)
                forward_edges[source].append(new_edge)
                reverse_edges[target].append(source)
                
                # Update sequences
                node_seqs[source] = new_seq
                
                changed = True
                break  # Restart after each merge as the graph changes

def write_merged_graph(filename, forward_edges, node_seqs):
    with open(filename, 'w') as f:
        for from_node, edges in forward_edges.items():
            for edge in edges:
                line = f"{from_node}\t{edge.to}\t{edge.weight}\t{edge.seq}\n"
                f.write(line)

def main():
    if len(sys.argv) != 3:
        print("Usage: python find_merges.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    forward_edges, reverse_edges, node_seqs = read_graph(input_file)
    
    print("Merging nodes...")
    merge_nodes(forward_edges, reverse_edges, node_seqs)
    
    print("Writing merged graph...")
    write_merged_graph(output_file, forward_edges, node_seqs)
    
    print(f"Merged graph written to {output_file}")

if __name__ == "__main__":
    main()