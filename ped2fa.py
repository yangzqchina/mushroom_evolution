# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 19:52:47 2016

@author: zqyang
"""
import sys

with open(sys.argv[2],'w') as fout:
    for line in open(sys.argv[1]):
        lines=line.rstrip().split('\t')
        head=lines[1]
        seq=''.join(lines[6:]).replace('0','N')
        fout.write(">%s\n%s\n"%(head,seq))
