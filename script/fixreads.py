#!/usr/bin/env python

import sys

fastareads=open(sys.argv[1],"r")
for line in fastareads:
    l=line[:-1]
    if len(l)>0 and l[0]==">":
#        sample = l[l.rfind(":")+1:]
        sample = sys.argv[2]
        if "barcodelabel=" in l:
            bpos = l.find(";barcodelabel")
            l = l[0:bpos]        
        print "%s;barcodelabel=%s" % (l,sample)
    else:
        print l
fastareads.close()
