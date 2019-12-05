#!/usr/bin/env python

import sys

ot = open(sys.argv[1],"r")
fa = open(sys.argv[2],"r")

seqIDs = set()
for line in fa:
    l = line.replace("\n","") 
    if l.startswith(">"):
        seqId = l[1:].split(" ")[0]
        seqIDs.add(seqId)  
    
fa.close()

firstLine = True
for line in ot:
    l = line.replace("\n","")
    if firstLine:
        firstLine = False
        print l
    else:
        otuID = l.split("\t")[0]
        if otuID in seqIDs:
            print l

