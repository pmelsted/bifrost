#!/usr/bin/python
#! -*- coding: utf-8 -*-

from sys import argv

if len(argv) < 3 or argv[2][-4:] != '.png':
    print "usage: ./make-graph.py <inputfile> <outputfile.png>"
    exit()


f = open(argv[1], "r").readlines()
d = {}
contigs = []
m = {}
numcontigs = len(f) / 2
for i in xrange(0, len(f), 2):
    info = f[i].replace("\n", "").strip().split(";")
    q = {}
    for item in info:
        if item[0] == '>':
            _id = int(item[7:])
        else:
            a, b = item.split(":")
            a, b = a.strip(), b.strip()
            if b:
                q[a] = map(int, b.split(","))
            else:
                q[a] = []
    c = f[i+1].replace("\n", "")
    q['contig'] = c
    q['ratio'] = q['Coveragesum'][0]  / ( 0.0 + q['Kmercount'][0])
    contigs.append(c)
    d[_id] = q
    m[_id] = c

it = d.items()
it.sort(key=lambda x: x[1]['ratio'], reverse=True)

import sys
sys.path.append('..')
sys.path.append('/usr/lib/graphviz/python/')
sys.path.append('/usr/lib64/graphviz/python/')
import gv

# Import pygraph
#from pygraph.classes.digraph import digraph
from pygraph.classes.graph import graph
from pygraph.readwrite.dot import write

# Graph creation
gr = graph()

#gr.add_nodes(m.values())

def p(_id):
    global d
    return "%d-%d-%d" % (_id, d[_id]['ratio'], d[_id]['Length'][0])

gr.add_nodes(map(p, range(numcontigs)))


for w in sorted(d.items()):
    k,v = w
    #print ">contig%d: %s" % (k, d[k]['contig'])
    for bw in v['Backwards']:
        if not gr.has_edge((p(bw), p(k))):
            gr.add_edge((p(bw), p(k)))
        #gr.add_edge((m[k], m[bw]))

    for fw in v['Forward']:
        if not gr.has_edge((p(k), p(fw))):
            gr.add_edge((p(k), p(fw)))
        #gr.add_edge((m[k], m[fw]))
#print d
# Draw as PNG
dot = write(gr)
gvv = gv.readstring(dot)
gv.layout(gvv,'dot')
gv.render(gvv,'png', argv[2])
