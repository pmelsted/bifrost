#!/usr/bin/python
#! -*- coding: utf-8 -*-
import sys

class Contig:
    def __init__(self, id, bases):
        self.id = id
        self.bases = bases

    def addinfo(self, length, ratio, bw, fw):
        self.length = length
        self.ratio = ratio
        self.bw = bw
        self.fw = fw


def createDict(prefix):
    non = lambda s : s.replace('\n', '').strip()
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

    for i in xrange(0, contigcount):
        assert clines[2*i] == ">contig%d" % i
        contigs.append(Contig(i, clines[1 + 2*i]))

    for i in xrange(contigcount):
        line = glines[1 + 3*i].split(" ")
        assert int(line[0]) == i
        length = int(line[1])
        ratio = float(line[2])
        bwcount = int(line[3])
        fwcount = int(line[4])
        bw = map(int, glines[2 + 3*i].split(" ")) if bwcount else []
        fw = map(int, glines[3 + 3*i].split(" ")) if fwcount else []
        contigs[i].addinfo(length, ratio, bw, fw)

    return contigs, KMERSIZE


def makeDot(contigs, KMERSIZE):
    s = ""
    s += "digraph G{\ngraph [rankdir=LR];\n node[shape=record]\n"
    max_cov = max(c.ratio for c in contigs)
    for c in contigs:
        now = c.bases
        if c.length >= 2*KMERSIZE:
            now = "%s .. (%d) .. %s" % (c.bases[:KMERSIZE], int(c.ratio), c.bases[-KMERSIZE:])
        s += '%s [style=filled, fillcolor=gray%s,label="%s"];\n'%(c.bases, int(100.0 - round(c.ratio / max_cov)), now)

    for c in contigs:
        for i in c.bw:
            s += "%s -> %s;\n" % (c.bases, contigs[i].bases)
        for i in c.fw:
            s += "%s -> %s;\n" % (c.bases, contigs[i].bases)

    s += "}"
    return s


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "usage: ./make-graph.py <prefix>"
        exit()

    prefix = sys.argv[1]
    contigs, KMERSIZE = createDict(prefix)
    dot = makeDot(contigs, KMERSIZE)
    out = open(prefix + ".dot", 'w')
    out.write(dot)
    out.close()

