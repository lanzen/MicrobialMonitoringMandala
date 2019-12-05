#!/usr/bin/env python

import sys

ot = open(sys.argv[2],"r")
fa = open(sys.argv[1],"r")
otuIDs = set()


firstLine = True
for line in ot:
    l = line.replace("\n","")
    if firstLine:
        firstLine = False
    else:
        otuID = l.split("\t")[0]
        otuIDs.add(otuID)

ot.close()
        
for line in fa:
    l = line.replace("\n","") 
    if l.startswith(">"):
        seqId = l[1:].split(" ")[0]
        if ";size" in seqId:
            seqId = seqId[:seqId.find(";size")]
        if seqId in otuIDs:
            otu_inc=True
            print ">%s" % seqId
        else:
            otu_inc=False
    elif otu_inc:
        print l
    
fa.close()

