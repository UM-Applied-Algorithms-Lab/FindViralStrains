from itertools import groupby
import argparse


def read_seq(file):
    with open(file, 'r') as f:
        lines = f.readlines()
##    ref = lines[0].split()
##    ref.pop(0) # remove name
##    primStrand = dict()
##    for i in ref:
##        n = i[:len(i)-1]
##        s = i[len(i)-1:]
##        primStrand[n] = (s == '+')
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
##        for k in range(0, len(nlist)):
##            pcount += 1 if (primStrand.get(nlist[k]) == (slist[k] == '+')) else 0
##            pcount -= 1 if (primStrand.get(nlist[k]) == (slist[k] != '+')) else 0
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

##def fasta_iter(fasta_name): # from https://www.biostars.org/p/710/
##    fh = open(fasta_name, 'r')
##    # ditch the boolean (x[0]) and just keep the header or sequence since
##    # we know they alternate.
##    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
##    for header in faiter:
##        # drop the ">"
##        headerStr = header.__next__()[1:].strip()
##        # join all sequence lines to one.
##        seq = "".join(s.strip() for s in faiter.__next__())
##        yield (headerStr, seq)

##def read_fa(file, edgeMers):
##    merCount = dict()
##    for m in edgeMers.values():
##        merCount[m] = 0
##    fiter = fasta_iter(file)
##    for ff in fiter:
##        headerStr, seq = ff
##        for m in edgeMers.values():
##            if m in seq or revcomp(m) in seq:
##                merCount[m] += 1
##    edgeCounts = []
##    for e,m in edgeMers.items():
##        edgeCounts.append([e, merCount[m]])
##    return edgeCounts

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
