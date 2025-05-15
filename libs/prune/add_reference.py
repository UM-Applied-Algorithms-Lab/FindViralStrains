import sys

dbg_file = sys.argv[1]
ref_file = sys.argv[2]
ref_and_dbg_file = sys.argv[3]

#Used in step 3
def reverse_complement(seq):
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

kmer_to_node = {}
seen_edges = set()

# Step 1: Parse existing graph
with open(dbg_file, 'r') as file, open(ref_and_dbg_file, 'w') as outfile:
    for line in file:
        parts = line.strip().split()
        if len(parts) != 4:
            print(f"Skipping malformed line: {line}")
            continue
        
        node1, node2, weight, edge_label = parts
        node1, node2 = int(node1), int(node2)        
        #If k changes from the existing 26
        k = len(edge_label) - 1
        src_kmer = edge_label[:k] #Equivalently [:-2]
        dst_kmer = edge_label[1:]

        kmer_to_node[src_kmer] = node1
        kmer_to_node[dst_kmer] = node2

        seen_edges.add((node1, node2))
        outfile.write(line)

# Step 2: Load the reference sequence
with open(ref_file, 'r') as file:
    lines = file.readlines()
    strain = "".join(line.strip() for line in lines[1:])

kmers = [strain[i:i+k] for i in range(len(strain) - k + 1)]

# Step 3: Assign new node IDs
node_id = [max(kmer_to_node.values(), default=0) + 1]

def get_or_assign_node_id(kmer):
    rc = reverse_complement(kmer)
    if kmer in kmer_to_node:
        return kmer_to_node[kmer]
    elif rc in kmer_to_node:
        return kmer_to_node[rc]
    else:
        kmer_to_node[kmer] = node_id[0]
        node_id[0] += 1
        return kmer_to_node[kmer]

with open(ref_and_dbg_file, 'a') as outfile:
    for i in range(len(kmers) - 1):
        src_kmer = kmers[i]
        dst_kmer = kmers[i + 1]
        src_id = get_or_assign_node_id(src_kmer)
        dst_id = get_or_assign_node_id(dst_kmer)

        edge_key = (src_id, dst_id)
        if edge_key not in seen_edges:
            seen_edges.add(edge_key)
            edge_label = src_kmer + dst_kmer[-1]
            outfile.write(f"{src_id}\t{dst_id}\t0\t{edge_label}\n")

