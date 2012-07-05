#!/usr/bin/python
#! -*- coding: utf-8 -*-
import sys

alpha = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
twin  = lambda x: ''.join(map(lambda z: alpha[z], list(x[::-1])))
non   = lambda s : s.replace('\n', '').strip()

class Contig:
    def __init__(self, id, bases):
        self.id = id
        self.bases = bases

    def addinfo(self, length, ratio, bw, fw):
        self.length = length
        self.ratio = ratio
        self.bw = bw
        self.fw = fw

    def __repr__(self):
        return "bw: %s fw: %s length: %d ratio %f seq: %s" % (self.bw, self.fw, self.length, self.ratio, self.bases)


def isNeighbour(seq1, seq2, KMERSIZE):
    aFirst = seq1[:KMERSIZE-1]
    bFirst = seq2[:KMERSIZE-1]
    aLast = seq1[-KMERSIZE+1:]
    bLast = seq2[-KMERSIZE+1:]
    if aLast == bFirst:
        return True
    elif aLast == twin(bLast):
        return True
    elif twin(aFirst) == bFirst:
        return True
    elif twin(aFirst) == twin(bLast):
        return True
    return False

def createDict(prefix):
    contigfile = prefix + ".contigs"
    graphfile = prefix + ".graph"
    try:
        clines = map(non, open(prefix + ".contigs", 'r').readlines())
    except:
        print "Could not find file %s.contigs, did you run ./BFGraph contigs -o %s ...?" % (prefix, prefix)
        exit()
    try:
        glines = map(non, open(prefix + ".graph", 'r').readlines())
    except:
        print "Could not find file %s.graph, did you run ./BFGraph contigs -o %s ...?" % (prefix, prefix)
        exit()
    print "Creating the De Brujin graph from the files %s and %s" % (contigfile, graphfile)

    contigcount, KMERSIZE = map(int, glines[0].split(" "))

    contigs = []

    d = {}
    for i in xrange(contigcount):
        assert clines[2*i] == ">contig%d" % i
        s = clines[1 + 2*i]

        first = s[:KMERSIZE-1]
        last = s[-KMERSIZE+1:]
        firsttwin = twin(first)
        lasttwin = twin(last)
        d.setdefault(first, []).append(i)
        d.setdefault(last, []).append(i)
        d.setdefault(firsttwin, []).append(i)
        d.setdefault(lasttwin, []).append(i)

        contigs.append(Contig(i, s))


    for i in xrange(contigcount):
        line = glines[1 + 3*i].split(" ")
        assert int(line[0]) == i
        length = int(line[1])
        ratio = float(line[2])
        bwcount = int(line[3])
        fwcount = int(line[4])
        bw = map(int, glines[2 + 3*i].split(" ")) if bwcount else []
        fw = map(int, glines[3 + 3*i].split(" ")) if fwcount else []
        c = contigs[i]
        c.addinfo(length, ratio, bw, fw)
        s = c.bases

        for back in bw:
            assert isNeighbour(contigs[back].bases, s, KMERSIZE)
        for forward in bw:
            assert isNeighbour(contigs[forward].bases, s, KMERSIZE)


    for c in contigs:
        s = c.bases
        first = s[:KMERSIZE-1]
        last = s[-KMERSIZE+1:]
        firsttwin = twin(first)
        lasttwin = twin(last)
        """
        print set(c.fw + c.bw + [c.id])
        print set(d[first] + d[last] + d[firsttwin] + d[lasttwin])
        """
        assert set(c.fw+c.bw+[c.id]) == set(d[first] + d[last] + d[firsttwin] + d[lasttwin])





    return contigs, KMERSIZE


def makeDot(contigs, KMERSIZE):
    s = "digraph G{\ngraph [rankdir=LR];\n node[shape=record]\n"
    max_cov = max(c.ratio for c in contigs)
    for c in contigs:
        now = c.bases
        if c.length >= 2*KMERSIZE:
            now = "%s .. (%d) .. %s" % (c.bases[:KMERSIZE], int(c.ratio), c.bases[-KMERSIZE:])
        s += '%s [style=filled, fillcolor=gray%s,label="%s"];\n'%(c.bases, int(100.0 - round(70*c.ratio / max_cov)), now)

    for c in contigs:
        for i in c.bw:
            s += "%s -> %s;\n" % (c.bases, contigs[i].bases)
        for i in c.fw:
            s += "%s -> %s;\n" % (c.bases, contigs[i].bases)

    s += "}"
    return s


def makeDot2(contigs, KMERSIZE):
    lines = ["digraph G {", "graph [rankdir=LR, fontcolor=red, fontname=\"Courier\"];", "node [shape=record];"]
    struct = {}
    score = {}
    i = 0
    for c in contigs:
        x = c.bases
        if x not in struct:
            struct[x] = (i,0)
            struct[twin(x)] = (i,1)
            i +=1
        score[x] = c.ratio
    max_score = max(score.values())
    for c in contigs:
        x = c.bases
        form = "style=filled, "
        form += "fillcolor=gray%s" % (int(100.0 - round(70*score[x]/max_score)),)
        lx = x
        ltx = twin(x)
        if c.length >= 2*KMERSIZE:
            lx = "%s .. (%d) .. %s" % (lx[:KMERSIZE], c.length, lx[-KMERSIZE:])
            ltx = "%s .. (%d) .. %s" % (ltx[:KMERSIZE], c.length, ltx[-KMERSIZE:])
        lines.append("%s[label=\"<%s> %s | <%s> %s\", %s];" % (struct[x][0],0,lx,1,ltx,form))


    done = {}
    for c in contigs:
        x = c.bases
        for fw in c.fw:
            o = contigs[fw].bases
            if (x,o) not in done:
                lines.append('%s:%s -> %s:%s;' % (struct[x][0],struct[x][1],struct[o][0],struct[o][1]))
                done[(x, o)] = 1
        for bw in c.bw:
            o = contigs[bw].bases
            if (o,x) not in done:
                lines.append('%s:%s -> %s:%s;' % (struct[o][0],struct[o][1],struct[x][0],struct[x][1]))
                done[(o, x)] = 1
    lines.append("}")
    return "\n".join(lines)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "usage: ./make-graph.py <prefix>"
        exit()

    prefix = sys.argv[1]
    contigs, KMERSIZE = createDict(prefix)
    dot = makeDot2(contigs, KMERSIZE)
    out = open(prefix + ".dot", 'w')
    out.write(dot)
    out.close()
