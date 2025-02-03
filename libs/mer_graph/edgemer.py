import argparse

def read_seq(file):
    """Reads the sequence file and yields edges one at a time."""
    with open(file, 'r') as f:
        for line in f:
            parts = line.split()
            parts.pop(0)  # remove name
            nlist = [j[:len(j)-1] for j in parts]
            slist = [j[len(j)-1:] for j in parts]
            pcount = 0
            for k in range(len(nlist) - 1):
                if pcount >= 0:  # vote is for primary strand
                    yield [nlist[k] + slist[k], nlist[k+1] + slist[k+1]]
                else:
                    yield [nlist[k+1] + ('+' if slist[k+1] == '-' else '-'), 
                           nlist[k] + ('+' if slist[k] == '-' else '-')]

def revcomp(seq):
    """Computes the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

def read_seg(file, edges, K):
    """Reads the segment file and yields edge-mer pairs one at a time."""
    nodeSeg = {}
    with open(file, 'r') as f:
        for line in f:
            l = line.split()
            nodeSeg[l[0] + '+'] = l[1]
            nodeSeg[l[0] + '-'] = revcomp(l[1])
    
    for edge in edges:
        mer = nodeSeg[edge[0]][-K] + nodeSeg[edge[1]][:K]
        yield edge[0] + ' ' + edge[1], mer

def write_graph(file, edges_func, edgeMers):
    """Writes the graph to the output file, processing data on-the-fly."""
    nodes = set()
    with open(file, 'w') as f:
        f.write('#graph 1\n')
        # First pass: count unique nodes
        for edge in edges_func():  # Call the generator function again
            nodes.add(edge[0])
            nodes.add(edge[1])
        f.write(str(len(nodes)) + '\n')
        # Second pass: write edge-mer pairs
        for e, m in edgeMers:
            f.write(e + ' ' + m + '\n')

def main():
    parser = argparse.ArgumentParser(prog='edgecnt',
                                     description='creates a count graph based on cuttlefish output and fasta')
    parser.add_argument('-k', '--kvalue', required=True, type=int, help='k-mer size')
    parser.add_argument('-c', '--cuttlefishprefix', required=True, help='prefix for cuttlefish output files')
    parser.add_argument('-o', '--output', required=True, help='output graph file')
    args = parser.parse_args()

    print("reading seq")
    edges_func = lambda: read_seq(args.cuttlefishprefix + '.cf_seq')  # Store generator function
    edges = edges_func()  # Create generator for the first iteration
    print("reading seg")
    edgeMers = read_seg(args.cuttlefishprefix + '.cf_seg', edges, args.kvalue)  # Generator for edge-mer pairs
    print("writing graph")
    write_graph(args.output, edges_func, edgeMers)  # Pass the generator function

if __name__ == '__main__':
    main()