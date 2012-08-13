#!/usr/bin/python
#! -*- coding: utf-8 -*-
import sys

alpha = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
twin  = lambda x : ''.join(map(lambda z: alpha[z], list(x[::-1])))

if __name__ == "__main__":
    if len(sys.argv) == 1:
        exit("Useage: python fix_contigs.py <input_file> <output_file>")
    inname = sys.argv[1]
    outname = sys.argv[2]
    f = open(inname, "r")
    cs = []
    for line in f.readlines():
        line = line.replace('\n', '')
        tmp = twin(line)
        if tmp < line:
            cs.append(tmp)
        else:
            cs.append(line)

    f = open(outname, "w")
    s = ""
    for line in sorted(cs):
        s += line+'\n'
    f.write(s)
    f.close()
