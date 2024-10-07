from itertools import groupby
import argparse


def read_seq(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    edges = []
    for i in lines:
        l = i.split()
        l.pop(0) # remove name
        nlist = []
        slist = []
        for j in l:
            nlist.append(j[:len(j)-1])
            slist.append(j[len(j)-1:])
        pcount = 0
        for k in range(0, len(nlist)-1):
            if pcount >= 0: # vote is for primary strand
                edge = [nlist[k] + slist[k], nlist[k+1] + slist[k+1]]
            else:
                edge = [nlist[k+1] + ('+' if slist[k+1] == '-' else '-'), nlist[k] + ('+' if slist[k] == '-' else '-')]
            edges.append(edge)
    return edges

def revcomp(seq):
    # complement strand
    rcseq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    rcseq = rcseq.upper()
    # reverse strand
    rcseq = rcseq[::-1]
    return rcseq

def read_seg(file, edges, K):
    nodeSeg = dict()
    with open(file, 'r') as f:
        for line in f:
            l = line.split()
            nodeSeg[l[0] + '+'] = l[1]
            nodeSeg[l[0] + '-'] = revcomp(l[1])
    edgeMers = dict()
    for edge in edges:
        mer = nodeSeg[edge[0]][len(nodeSeg[edge[0]]) - K] + nodeSeg[edge[1]][:K]
        edgeMers[edge[0] + ' ' + edge[1]] = mer
    return edgeMers

def write_graph(file, edges, edgeMers):
    nodes = set()
    for e in edges: # only consider nodes that belong to edges
        nodes.add(e[0])
        nodes.add(e[1])
    with open(file, 'w') as f:
        f.write('#graph 1\n')
        f.write(str(len(nodes)) + '\n')
        for e,m in edgeMers.items():
            # f.write(e + ' ' + '1' + '\n')
            f.write(e + ' ' + m + '\n')

def main():
    parser = argparse.ArgumentParser(prog='edgecnt',
                                     description='creates a count graph based on cuttlefish output and fasta')
    parser.add_argument('-k', '--kvalue', required=True)
    parser.add_argument('-c', '--cuttlefishprefix', required=True)
    #parser.add_argument('-f', '--fasta', required=True)
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    print("reading seq")
    edges = read_seq(args.cuttlefishprefix + '.cf_seq')
    print("reading seg")
    edgeMers = read_seg(args.cuttlefishprefix + '.cf_seg', edges, int(args.kvalue))
    print("reading fasta, counting edgemers")
    #edgeCounts = read_fa(args.fasta, edgeMers)
    print("writing graph")
    write_graph(args.output, edges, edgeMers)

if __name__ == '__main__':
    main()
