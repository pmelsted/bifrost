#!/usr/bin/python
#! -*- coding: utf-8 -*-
import sys
try:
    import pydot
except Exception as e:
    print "\nERROR: You do not have pydot installed, try: sudo pip install pydot\n"
    raise(e)

def writeToPNG(dotfile, pngfile):
    try:
        f = open(dotfile, 'r')
    except:
        print "Could not find file %s" % dotfile
        exit()

    try:
        g = pydot.graph_from_dot_data(f.read())
    except:
        print "\nERROR: Are you sure you have graphviz installed? You could try: sudo apt-get install graphviz\n"
        raise
    print "Writing the graph to %s" % pngfile
    g.write(pngfile, format='png')

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "usage: ./dot2png.py <file>"
        exit()
    if sys.argv[1][-4:] == ".dot":
        infile = sys.argv[1]
        prefix = sys.argv[1][0:-4]
    else:
        prefix = sys.argv[1]
        infile = prefix + ".dot"

    outfile = prefix + ".png"
    writeToPNG(infile, outfile)
