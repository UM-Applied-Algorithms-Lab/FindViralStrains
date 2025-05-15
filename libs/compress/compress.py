import sys
import re
from collections import defaultdict

def read_graph(filename):
    edges = defaultdict(list)
    node_seqs = {}
    with open(filename, 'r') as f:
        for line in f:
            # Skip empty lines or comment lines
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Split on any whitespace (tabs or multiple spaces)
            parts = re.split(r'\s+', line)
            if len(parts) < 4:
                continue
                
            try:
                from_node = int(parts[0])
                to_node = int(parts[1])
                weight = int(parts[2])
                seq = ' '.join(parts[3:])  # Join remaining parts
                edges[from_node].append((to_node, weight, seq))
                node_seqs[from_node] = seq
                if to_node not in node_seqs:
                    node_seqs[to_node] = ""
            except (ValueError, IndexError) as e:
                print(f"Warning: Skipping malformed line: {line} (Error: {e})")
                continue
    return edges, node_seqs

def compress_unitigs(edges, node_seqs, exclude_nodes={0, 1}):
    visited = set()
    unitigs = []
    total_compressed = 0

    for node in list(edges.keys()):
        if node in visited or node in exclude_nodes:
            continue

        # Start new unitig
        current_node = node
        unitig_nodes = [current_node]
        unitig_seqs = [node_seqs.get(current_node, "")]
        visited.add(current_node)

        # Extend forward
        while True:
            neighbors = edges.get(current_node, [])
            if len(neighbors) != 1:
                break
            next_node, weight, seq = neighbors[0]
            if next_node in exclude_nodes or next_node in visited:
                break
            # Check incoming edges
            incoming = sum(1 for n in edges if any(to == next_node for to, _, _ in edges[n]))
            if incoming != 1:
                break
            unitig_nodes.append(next_node)
            unitig_seqs.append(seq)
            visited.add(next_node)
            current_node = next_node

        # Combine sequences
        if len(unitig_nodes) == 1:
            unitig_seq = node_seqs.get(unitig_nodes[0], "")
        else:
            unitig_seq = unitig_seqs[0]
            for seq in unitig_seqs[1:]:
                if seq:
                    unitig_seq += seq[-1]

        if len(unitig_nodes) > 1:
            total_compressed += len(unitig_nodes) - 1

        unitigs.append((unitig_nodes[0], unitig_nodes[-1], len(unitig_nodes), unitig_seq))

    # Handle isolated nodes
    all_nodes = set(edges.keys())
    all_nodes.update(to for neighbors in edges.values() for to, _, _ in neighbors)
    for node in all_nodes:
        if node not in visited and node not in exclude_nodes:
            unitigs.append((node, node, 1, node_seqs.get(node, "")))
        elif node in exclude_nodes and node in node_seqs:
            unitigs.append((node, node, 1, node_seqs[node]))

    return unitigs, total_compressed

def write_unitigs(unitigs, output_filename, total_compressed):
    with open(output_filename, 'w') as f:
        f.write(f"# Total nodes compressed: {total_compressed}\n")
        for start, end, length, seq in sorted(unitigs, key=lambda x: x[0]):
            f.write(f"{start}\t{end}\t{length}\t{seq}\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python compress_unitigs.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    try:
        print(f"Processing {input_file}...")
        edges, node_seqs = read_graph(input_file)
        print(f"Graph loaded with {len(edges)} nodes")
        unitigs, total_compressed = compress_unitigs(edges, node_seqs)
        write_unitigs(unitigs, output_file, total_compressed)
        print(f"Success! Compressed {total_compressed} nodes into {len(unitigs)} unitigs.")
        print(f"Results written to {output_file}")
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()