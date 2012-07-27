#!/usr/bin/python
#! -*- coding: utf-8 -*-
import sys

alpha = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
twin  = lambda x : ''.join(map(lambda z: alpha[z], list(x[::-1])))
rep   = lambda s : min(s, twin(s))
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
        return "id: %d bw: %s fw: %s length: %d ratio %f seq: %s" % (self.id, self.bw, self.fw, self.length, self.ratio, self.bases)


def isNeighbour(seq1, seq2, KMERSIZE):
    aFirst = seq1[:KMERSIZE-1]
    bFirst = seq2[:KMERSIZE-1]
    aLast = seq1[-KMERSIZE+1:]
    bLast = seq2[-KMERSIZE+1:]
    if aLast == bFirst:
        return (0,0)
    elif aLast == twin(bLast):
        return (0,1)
    elif twin(aFirst) == bFirst:
        return (1,0)
    elif twin(aFirst) == twin(bLast):
        return (1,1)
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

    firsts = {}
    lasts = {}
    for i in xrange(contigcount):
        assert clines[2*i] == ">contig%d" % i
        s = clines[1 + 2*i]
        contigs.append(Contig(i, s))

        # For asserting that the fw and bw contigs are right
        first = s[:KMERSIZE]
        last = s[-KMERSIZE:]
        firsts.setdefault(first, []).append(i)
        lasts.setdefault(last, []).append(i)



    for i in xrange(contigcount):
        line = glines[1 + 5*i].split("_")
        assert int(line[0]) == i
        length = int(line[1])
        ratio = float(line[2])
        bw = map(int, glines[2 + 5*i].split(" ")) if glines[2 + 5*i] else []
        fw = map(int, glines[3 + 5*i].split(" ")) if glines[3 + 5*i] else []
        ibw = map(int, glines[4 + 5*i].split(" ")) if glines[4+ 5*i] else []
        ifw = map(int, glines[5 + 5*i].split(" ")) if glines[5 + 5*i]  else []
        c = contigs[i]
        c.addinfo(length, ratio, bw, fw)
        s = c.bases

        for back in bw:
            if not isNeighbour(contigs[back].bases, s, KMERSIZE):
                print i, s, bw, "Backward to %d" % back
                assert False
        for forward in fw:
            if not isNeighbour(contigs[forward].bases, s, KMERSIZE):
                print i, s, fw, "Forward to %d" % forward
                assert False


    for c in contigs:
        s = c.bases
        first = s[:KMERSIZE]
        last = s[-KMERSIZE:]


        #print ">contig%d\n%s" % (c.id, s)

        _fws = set()
        _bws = set()
        for char in "ACGT":
            fwkm = last[1:] + char
            bwkm = char + first[:-1]
            _fws.update(firsts.get(fwkm, []))
            _fws.update(lasts.get(twin(fwkm), []))
            _bws.update(lasts.get(bwkm, []))
            _bws.update(firsts.get(twin(bwkm), []))
        if not _fws == set(c.fw):
            print c, _fws, "Forward"
            assert False

        if not _bws == set(c.bw):
            print c, _bws, "Backward"
            assert False

    return contigs, KMERSIZE


def makeDot(contigs, KMERSIZE):
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
        for nbr in c.fw + c.bw:
            o = contigs[nbr].bases
            if (x,o) not in done:
                a, b = isNeighbour(x, o, KMERSIZE)
                lines.append('%s:%s -> %s:%s;' % (struct[x][0],a,struct[o][0], b))
                done[(x, o)] = 1
                done[(o, x)] = 1
    lines.append("}")
    return "\n".join(lines)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "usage: ./make-graph.py <prefix>"
        exit()

    prefix = sys.argv[1]
    contigs, KMERSIZE = createDict(prefix)
    dot = makeDot(contigs, KMERSIZE)
    out = open(prefix + ".dot", 'w')
    out.write(dot)
    print "Wrote the dot graph to %s" % (prefix + ".dot")
    out.close()
